function [varargout] = plot_contour(xy, z_data, amp, new_fig, flipped, is_cont, varargin)
%% Description
% Creates a contour plot from a set of 2D input locations (mesh).
% Automatically creates a grid in order for the contour plot to work.
% Inputs:
% - mesh: matrix of size Mx2. Rows are observations, columns are x- and
% y-locations.
% - y_data: column vector of size M. Data to plot.

%% Input
if nargin < 2
    error('You must provide two input arguments.');
end
if nargin < 3 || isempty(amp)
    amp = 1;
end
if nargin < 4 || isempty(new_fig)
    new_fig = true;
end
if nargin < 5 || isempty(flipped)
    flipped = false;
end
if nargin < 6 || isempty(is_cont)
    is_cont = false;
end
%% Pre-processing
% Get unique coordinates
x = unique(sort(xy(:,1))); 
y = unique(sort(xy(:,2))); 
% Make mesh denser
if amp > 1
    % Make the mesh more dense 
    deltax = abs(x(2) - x(1));
    deltay = abs(y(2) - y(1));
    delta = min([deltax, deltay]) / amp;
    x = x(1) - delta : delta : x(end) + delta;
    y = y(1) - delta : delta : y(end) + delta;
end
x = x(:);
y = y(:);
% Create variables to plot
method = 'nearest';
[xx, yy] = meshgrid(x, y);
F = scatteredInterpolant(xy(:,1), xy(:,2), z_data(:), method);
zz = F(xx, yy);
if flipped
    F = scatteredInterpolant(xy(:,1), xy(:,2), z_data, method);
    zz2 = F(xx, yy);
    zz = [flip(zz,2), zz2];
    [xx, yy] = meshgrid([flip(-x); x], y);
end
%% Plot
if new_fig
    figure();
end
if is_cont
    contourf(xx, yy, zz);
else
    if flipped
        imshow(zz, [], 'XData', x, 'YData', y, 'Colormap', jet);
        axis xy; % Flip the figure
        axis on;
    else
        imshow(zz, [], 'XData', x, 'YData', y, 'Colormap', jet);
        axis xy; % Flip the figure
        axis on; 
    end
end
xlabel('x [m]');
ylabel('y [m]');
end



% Make the mesh more dense 
% deltax = abs(x(2) - x(1));
% deltay = abs(y(2) - y(1));
% delta = min([deltax, deltay]) / 3.0;
% x = min((obj.mesh(:,1))) : delta : max((obj.mesh(:,1))); x = x(:);
% y = min((obj.mesh(:,2))) : delta : max((obj.mesh(:,2))); y = y(:);
% % Get data to plot
% [xx, yy] = meshgrid(x, y);
% F = scatteredInterpolant(obj.mesh(:,1), obj.mesh(:,2), y2plot, method);
% zz = F(xx, yy);
% F = scatteredInterpolant(obj.mesh(:,1), obj.mesh(:,2), y_data, method);
% zz2 = F(xx, yy);
% zz = [flip(zz,2), zz2];
% [xx, yy] = meshgrid([flip(-x); x], y);
% % Plot
% figure();
% if cont
%     y_min = min(y2plot); y_max = max(y2plot); % Range
%     v = y_min + (0 : .05 : 1) * (y_max - y_min); % Color values
%     [~, h, ~] = contourf(xx, yy, zz, v); 
%     set(h(:), 'LineStyle', 'none');
% else
%     imshow(zz, [], 'XData', [flip(-x); x], 'YData', y, 'Colormap', jet);
%     axis xy; % Flip the figure
%     axis on; 
% end
