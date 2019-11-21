function [varargout] = plot_reconstructionErrors(obj, varargin)
% DESCRIPTION
% varargin{1}: model
% varargin{2}: var
% varargin{3}: plot_simple
% 

% Check the points are there 
if isempty(obj.training_points) || isempty(obj.prediction_points)
    fprintf('No training or prediction points are available.\n');
    cond_no_points = true;
else
    cond_no_points = false;
end

% Check the dimension of the input space
if obj.is_mesh_variable
    x_train = obj.training_points;
    x_test = obj.prediction_points;
    dim = size(x_train, 2);
else
    x_train = unique( obj.training_points(:,size(obj.mesh,2)+1:end) );
    x_test = unique( obj.prediction_points(:,size(obj.mesh,2)+1:end) );
    dim = size(obj.training_points, 2) - size(obj.mesh,2);
end
if dim > 3 
    fprintf('Dimension of input space is > 3.\n');
    cond_no_points = true;
end

% Input
n_args = length(varargin);
model = 'pca';
if n_args > 0
    model = varargin{1};
    I = strcmp({'pca', 'lpca', 'cpca', 'lcpca', 'pca_var', 'lpca_var', 'cpca_var', 'lcpca_var'}, model);
    if sum(I) < 1
        error('Model not valid.')
    end
end

% Get correct errors
if strcmp(model, 'pca')
    y_err = obj.pca_reconstruction_error_observations;
elseif strcmp(model, 'lpca')
    y_err = obj.local_pca_reconstruction_error_observations;
elseif strcmp(model, 'cpca')
    y_err = obj.cpca_reconstruction_error_observations;
elseif strcmp(model, 'lcpca')
    y_err = obj.local_cpca_reconstruction_error_variables;
elseif strcmp(model, 'pca_var')
    y_err = obj.pca_reconstruction_error_variables;
elseif strcmp(model, 'lpca_var')
    y_err = obj.local_pca_reconstruction_error_variables;
elseif strcmp(model, 'cpca_var')
    y_err = obj.cpca_reconstruction_error_variables;
elseif strcmp(model, 'lcpca_var')
    y_err = obj.local_cpca_reconstruction_error_variables;
end

% If present, run only for one specific variable
var = 'all';
if n_args > 1 && ~isempty(varargin{2})
    var = varargin{2};
end
if ~strcmp(var, 'all')
    y_err = obj.get_variable_errors(false, var, model, obj.is_mesh_variable);
end

% Number of points 
m = size(x_train, 1);

% Symbols
s1 = 'ko';
s2 = 'bx';

% Plot variables?
is_var = strcmp(model(end-2:end), 'var');

% Choose plotting method
plot_simple = false;
if n_args > 2
    plot_simple = varargin{3};
end

% Start new figure
figure();
if ~is_var && ~cond_no_points && ~plot_simple
    % Input space with dimension = 1
    if dim == 1
        plot(x_train, 100 * y_err, 'LineStyle','--', 'Marker','*', 'Color','k');  
        grid on; hold on;
        ylabel('$\mathrm{Error \ [\%]}$','Interpreter','LaTex');
    end

    % Input space with dimension = 2
    if dim == 2
        plot(x_train(:,1), x_train(:,2), 'ko');
        hold on; grid on;
        plot(x_test(:,1), x_test(:,2), s2);
        hold on; grid on;
        % Plot errors
        for ii = 1 : m
            x = x_train(ii,1);
            y = x_train(ii,2);
            txt = [' ', num2str(100 * y_err(ii),3), '%'];
            text(x,y,txt);
        end
    end

    % Input space with dimension = 3
    if dim == 3 
        plot3(x_train(:,1), x_train(:,2), x_train(:,3), 'ko');
        hold on; grid on;
        plot3(x_test(:,1), x_test(:,2), x_test(:,3), s2);
        hold on; grid on;
        % Plot errors
        for ii = 1 : m
            x = x_train(ii,1);
            y = x_train(ii,2);
            z = x_train(ii,3);
            txt = [' ', num2str(100 * y_err(ii),3), '%'];
            text(x,y,z,txt);
        end
    end
else
    plot(1:length(y_err), 100 * y_err, 'LineStyle','--', 'Marker','*', 'Color','k');  
    grid on; hold on;
    xlabel('Input point index $[-]$','Interpreter','LaTex');
    ylabel('Error $[\%]$','Interpreter','LaTex');
end

% Legend and title
if strcmp(var, 'all')
    tit = ['\textbf{Reconstruction errors: ', model, '}'];
else
    tit = ['\textbf{Reconstruction errors (', var, '): ', model, '}'];
end
temp = obj.pca_approximation_order;
tit = [tit, '; ', num2str(temp), ' PCs']; 
if strcmp(model(1), 'l')
    tit = [tit, ' ', num2str(obj.number_of_clusters)];
end
title(tit,'Interpreter','LaTex');
hold off;

end

