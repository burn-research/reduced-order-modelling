function [varargout] = plot_Kriging_surface(training_points, prediction_points, values, kriged_values, varargin)
%% Description:
% Plot a genral Kriging response surface.

%% Input
% Label
label_string = [];
if ~isempty(varargin)
    label_string = varargin{1};
end
% Knwon targets
known_targets = [];
if length(varargin) > 1
    known_targets = varargin{2};
end

dim = size(training_points, 2);

%% Plot
% Plot routines
if dim == 1
    plot1D(training_points, prediction_points, values, kriged_values, label_string, known_targets);
elseif dim == 2
    plot2D(training_points, prediction_points, values, kriged_values, label_string, known_targets);
elseif dim == 3
    plot3D(training_points, prediction_points, values, kriged_values, label_string, known_targets);
else
    error('Not possible to plot.');
end

if nargout > 0
    varargout{1} = h;
end

end


function [varargout] = plot1D(training_points, prediction_points, values, kriged_values, label_string, varargin)

% Knwon targets
known_targets = [];
if ~isempty(varargin)
    known_targets = varargin{1};
end

% Plot
plot(training_points, values, '*b');
hold on;
plot(prediction_points, kriged_values, 'or');
if ~isempty(known_targets)
    plot(prediction_points, kriged_values, '+');
end
grid on;
hold off;
title('Response Surface');
if ~isempty(label_string)
    ylabel(label_string);
end

if nargout > 0
    varargout{1} = h;
end

end 


function [varargout] = plot2D(training_points, prediction_points, values, kriged_values, label_string, varargin)

% Knwon targets
known_targets = [];
if ~isempty(varargin)
    known_targets = varargin{1};
end

if ~isempty(kriged_values)
    [points, idx_points] = sortrows([training_points; prediction_points]);
else
    [points, idx_points] = sortrows(training_points);
end

% Domain
x = unique(points(:,1));
y = unique(points(:,2));

% Make the mesh more dense 
dx = abs(x(2) - x(1));
dy = abs(y(2) - y(1));
fac = 1 / min([length(x), length(y)]);
delta = min([dx, dy]) * fac;
x = points(1,1) : delta : points(end,1);
y = points(1,2) : delta : points(end,2);

% Meshgrid
[xx, yy] = meshgrid(x, y);

% z-axis values
if ~isempty(kriged_values)
    y_data = [values, kriged_values]';
else
    y_data = values';
end
y_data = y_data(idx_points);

% Plot
if strcmp(computer, 'MACI64')
    F = scatteredInterpolant(points(:,1), points(:,2), y_data, 'natural');
    zz = F(xx, yy);
    surf(xx, yy, zz, 'linestyle', 'none'); hold on;
    alpha(0.85);
end
plot3(training_points(:,1), training_points(:,2), values, '*b'); hold on;
plot3(prediction_points(:,1), prediction_points(:,2), kriged_values, 'or'); hold on;
if ~isempty(known_targets)
    plot3(prediction_points(:,1), prediction_points(:,2), known_targets, 'x'); hold on;
end
grid on;
hold off;
title('Response Surface');
if ~isempty(label_string)
    zlabel(label_string);
end

if nargout > 0
    varargout{1} = h;
end

end 


function [varargout] = plot3D(training_points, prediction_points, values, kriged_values, label_string, varargin)

scatter3(training_points(:,1), training_points(:,2), training_points(:,3), 44, ...
    values, 'filled');
hold on;
scatter3(prediction_points(:,1), prediction_points(:,2), prediction_points(:,3), 44, ...
    kriged_values, 'filled');
grid on;
hold off;
title('Response Surface');
if ~isempty(label_string)
    zlabel(label_string);
end

end






