function [varargout] = plot_contour(obj, variable, point, varargin)
% DESCRIPTION: Works only for Zhy Li's case.

if ~iscell(obj.mesh) && size(obj.mesh,2) ~= 2
    error('It looks like there is no 2D mesh.');
elseif iscell(obj.mesh) && size(obj.mesh{1},2) ~= 2
    error('It looks like there is no 2D mesh.');
end

if size(point,1) > 1
    error('Too many points provided.');
end

if ~isempty(varargin)
    model = varargin{1};
else
    model = 'data';
end

% Get vectors to plot
if strcmp(model, 'pc')
    y_data = obj.get_variable(point, variable, 'pc');
else
    y_data = obj.get_variable(point, variable, 'data');
end
if ischar(variable)
    y2plot = obj.get_variable(point, variable, model);
else
    % Possibility to input the data to plot
    y2plot = variable;
end

method = 'natural';
cont = false;

% Get unique coordinates
if iscell(obj.mesh)
    % Find variable's name's position
    for ii = 1 : length(obj.variable_names)
        if strcmp(variable, obj.variable_names{ii})
            var_i = ii;
        end
    end
    x = unique(obj.mesh{var_i}(:,1));
    y = unique(obj.mesh{var_i}(:,2));
else
    x = unique(obj.mesh(:,1));
    y = unique(obj.mesh(:,2)); 
end

if size(obj.mesh, 1) > 13130
    try
        % Make the mesh more dense 
        idx_half = round(size(obj.mesh,1)/2); % Take half mesh
        deltax = abs(x(2) - x(1));
        deltay = abs(y(2) - y(1));
        delta = min([deltax, deltay]) / 6.0;
        x = min((obj.mesh(1:idx_half,1))) : delta : max((obj.mesh(1:idx_half,1)));
        y = min((obj.mesh(:,2))) : delta : max((obj.mesh(:,2)));
        x2 = min((obj.mesh(1+idx_half:end,1))) : delta : max((obj.mesh(1+idx_half:end,1)));
        % Get data to plot
        [xx, yy] = meshgrid(x, y);
        F = scatteredInterpolant(obj.mesh(:,1), obj.mesh(:,2), y2plot, method);
        zz = F(xx, yy);
        [xx2, yy2] = meshgrid(x2, y);
        F = scatteredInterpolant(obj.mesh(:,1), obj.mesh(:,2), y_data, method);
        zz2 = F(xx2, yy2);
        [xx, yy] = meshgrid([x(:); x2(:)], y);
        zz = [zz, zz2];
        % Plot
        figure();
        if cont
            y_min = min(y2plot); y_max = max(y2plot); % Range
            v = y_min + (0 : .05 : 1) * (y_max - y_min); % Color values
            [~, h, ~] = contourf(xx, yy, zz, v); 
            set(h(:), 'LineStyle', 'none');
        else
            imshow(zz, [], 'XData', [x(:); x2(:)], 'YData', y, 'Colormap', jet);
            axis xy; % Flip the figure
            axis on; 
        end
    catch ME
        figure(); grid on;
        scatter(obj.mesh(:,1), obj.mesh(:,2), 35, y2plot, 'filled');
    end
else
    % Make the mesh more dense 
    deltax = abs(x(2) - x(1));
    deltay = abs(y(2) - y(1));
    delta = min([deltax, deltay]) / 1;
    x = min((obj.mesh(:,1))) : delta : max((obj.mesh(:,1))); x = x(:);
    y = min((obj.mesh(:,2))) : delta : max((obj.mesh(:,2))); y = y(:);
    % Get data to plot
    [xx, yy] = meshgrid(x, y);
    F = scatteredInterpolant(obj.mesh(:,1), obj.mesh(:,2), y2plot, method);
    zz = F(xx, yy);
    F = scatteredInterpolant(obj.mesh(:,1), obj.mesh(:,2), y_data, method);
    zz2 = F(xx, yy);
    zz = [flip(zz,2), zz2];
    [xx, yy] = meshgrid([flip(-x); x], y);
    % Plot
    figure();
    if cont
        y_min = min(y2plot); y_max = max(y2plot); % Range
        v = y_min + (0 : .05 : 1) * (y_max - y_min); % Color values
        [~, h, ~] = contourf(xx, yy, zz, v); 
        set(h(:), 'LineStyle', 'none');
    else
        imshow(zz, [], 'XData', [flip(-x); x], 'YData', y, 'Colormap', jet);
        axis xy; % Flip the figure
        axis on; 
    end
    set(gca, 'XDir','reverse');
end
% scatter(obj.mesh(:,1), obj.mesh(:,2), 15, y2plot, 'filled');
xlabel('x [m]');
ylabel('y [m]');
title([variable, ': data - ', model]);
set(findall(gcf,'-property','FontSize'),'FontSize',25);
if nargout > 0
    varargout{1} = zz;
end

end



% Local functions
function idx = get_var_dom(A, n, var_i)

idx = A(1+n*(var_i-1):n*var_i, :);

end






