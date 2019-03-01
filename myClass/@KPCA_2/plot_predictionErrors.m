function [varargout] = plot_predictionErrors(obj, varargin)

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
    I = strcmp({'pca', 'lpca', 'cpca', 'lcpca', 'direct', 'pca_var', 'lpca_var', 'cpca_var', 'lcpca_var', 'direct_var'}, model);
    if sum(I) < 1
        error('Model not valid.')
    end
end

% Get correct errors
if strcmp(model, 'pca')
    y_err = obj.kpca_prediction_error_observations;
elseif strcmp(model, 'lpca')
    y_err = obj.klpca_prediction_error_observations;
elseif strcmp(model, 'cpca')
    y_err = obj.kcpca_prediction_error_observations;
elseif strcmp(model, 'lcpca')
    y_err = obj.klcpca_prediction_error_observations;
elseif strcmp(model, 'direct')
    y_err = obj.direct_prediction_error_observations;
elseif strcmp(model, 'pca_var')
    y_err = obj.kpca_prediction_error_variables;
elseif strcmp(model, 'lpca_var')
    y_err = obj.klpca_prediction_error_variables;
elseif strcmp(model, 'cpca_var')
    y_err = obj.kcpca_prediction_error_variables;
elseif strcmp(model, 'lcpca_var')
    y_err = obj.klcpca_prediction_error_variables;
elseif strcmp(model, 'direct_var')
    y_err = obj.direct_prediction_error_variables;
end

% If present, run only for one specific variable
var = 'all';
if n_args > 1
    var = varargin{2};
end
if ~strcmp(var, 'all')
    y_err = obj.get_variable_errors(true, var, model, obj.is_mesh_variable);
end

% Number of points 
m = size(x_test, 1);

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
%         [x, I] = sortrows([x_train; x_test]);
%         y = [0*x_train; y_err]; y = y(I);
        x = x_test;
        y = y_err;
        plot(x, 100 * y, 'LineStyle','--', 'Marker','*', 'Color','k'); 
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
            x = x_test(ii,1);
            y = x_test(ii,2);
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
            x = x_test(ii,1);
            y = x_test(ii,2);
            z = x_test(ii,3);
            txt = [' ', num2str(100 * y_err(ii),3), '%'];
            text(x,y,z,txt);
        end
    end
else
    plot(1:length(y_err), 100 * y_err, 'LineStyle','--', 'Marker','*', 'Color','k');  
    grid on; hold on;
    xlabel('$\mathrm{Input \ point \ index \ [-]}$','Interpreter','LaTex');
    ylabel('$\mathrm{Error \ [\%]}$','Interpreter','LaTex');
end

% Legend and title
k = obj.pca_approximation_order;
if strcmp(var, 'all')
    tit = ['\textbf{Prediction errors: ', model, ' ', num2str(k), ' PCs}'];
else
    tit = ['\textbf{Prediction errors (', var, '): ', model, ' ', num2str(k), ' PCs}'];
end
if strcmp(model(1), 'l')
    tit = [tit, '; Clusters: ', num2str(obj.number_of_clusters)];
end
title(tit,'Interpreter','LaTex');
hold off;

end


