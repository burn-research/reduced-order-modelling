function [varargout] = plot_inputSpaceErrors(obj, varargin)


% Check the points are there 
if isempty(obj.training_points) || isempty(obj.prediction_points)
    fprintf('No training or prediction points are available.\n');
    return;
end

% Check the dimension of the input space
dim = size(obj.training_points, 2);
if dim > 3
    fprintf('Dimension of input space is > 3.\nNo plot possible.\n');
    return;    
end

% Input
model = 'pca';
if ~isempty(varargin)
    model = varargin{1};
end

% Model
if strcmp(model, 'pca')
    z_xy = obj.kpca_prediction_error_observations;
elseif strcmp(model, 'lpca')
    z_xy = obj.klpca_prediction_error_observations;
elseif strcmp(model, 'cpca')
    z_xy = obj.kcpca_prediction_error_observations;
elseif strcmp(model, 'lcpca')
    z_xy = obj.klcpca_prediction_error_observations;
end

% Number of points 
m = size(obj.prediction_points, 1);

% Symbols
s1 = 'ko'; % Marker for training points
s2 = 'bx'; % Marker for prediction points

% Start new figure
figure();

% Input space with dimension = 1
if dim == 1
    y_vec = obj.training_points * 0 + 1;
    plot(obj.training_points, y_vec, s1);
    hold on; grid on;
    plot(obj.prediction_points, y_vec, s2);
    hold on; grid on;
    % Plot errors
    for ii = 1 : m
        x = obj.prediction_points(ii);
        y = 1;
        txt = [' ', num2str(100 * z_xy(ii),3), '%'];
        text(x,y,txt);
    end
end

% Input space with dimension = 2
if dim == 2
    plot(obj.training_points(:,1), obj.training_points(:,2), s1);
    hold on; grid on;
    plot(obj.prediction_points(:,1), obj.prediction_points(:,2), s2);
    hold on; grid on;
    % Plot errors
    for ii = 1 : m
        x = obj.prediction_points(ii,1);
        y = obj.prediction_points(ii,2);
        txt = [' ', num2str(100 * z_xy(ii),3), '%'];
        text(x,y,txt);
    end
end

% Input space with dimension = 3
if dim == 3 
    plot3(obj.training_points(:,1), obj.training_points(:,2), obj.training_points(:,3), s1);
    hold on; grid on;
    plot3(obj.prediction_points(:,1), obj.prediction_points(:,2), obj.prediction_points(:,3), s2);
    hold on; grid on;
    % Plot errors
    for ii = 1 : m
        x = obj.prediction_points(ii,1);
        y = obj.prediction_points(ii,2);
        z = obj.prediction_points(ii,3);
        txt = [' ', num2str(100 * z_xy(ii),3), '%'];
        text(x,y,z,txt);
    end
end

% Legend and title
title(['Errors for prediction points: ', model]);
hold off;

end


