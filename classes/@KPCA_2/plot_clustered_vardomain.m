function [varargout] = plot_clustered_vardomain(obj, varargin)
%% Description:
% Plot how the mesh-points have been clustered for a specific variable.
% Input
    % (1) Physical variable to be analysed, or input space point.
%

%% Input
% Mesh size
[n_points, m_dim] = size(obj.mesh);

% Check this condition
if m_dim > 2
    warning('Not possible to plot this with a mesh whose dimension is higher than 2.');
    return
end

% We need info on the input parameters
if ~obj.is_mesh_variable
    off_set = m_dim;
else
    off_set = 0;
end
p_dim = size(obj.training_points,2) - off_set; % Get dimension
p = cell(p_dim, 1);
p_all = obj.training_points(:, off_set+1:end);

n_args = length(varargin);
% Input (1): 
if n_args > 0
    % User-input
    if obj.is_mesh_variable
        if ~isa(varargin{1}, 'char')
            arg_1 = obj.variable_names{varargin{1}};
        else
            arg_1 = varargin{1};
        end
    else
        if ischar(varargin{1})
            error('You were probably trying to select a variable but mesh is a parameter.');
        end
        arg_1 = varargin{1};
    end
else
    % Default
    if obj.is_mesh_variable
        arg_1 = obj.variable_names{1};
    else
        arg_1 = p_all(1,:);
    end
end 

% Input (2):
if n_args > 1
    var = varargin{2};
else
    var = obj.variable_names{1};
end


%% Main
% Find variable's name's position (works on arg_1)
if ischar(arg_1)
    for ii = 1 : length(obj.variable_names)
        if strcmp(arg_1, obj.variable_names{ii})
            var_i = ii;
        end
    end
else
    n_obs = round(size(p_all,1) / n_points); % Number of observations
    d = ones(n_obs, 1);
    val = zeros(n_obs, 1);
    for i = 1 : n_obs
        val(i) = p_all(n_points*(i-1)+1,:);
        temp = (val(i) - arg_1).^2;
        d(i) = sqrt(sum(temp, 2));
    end
    [~, var_i] = min(d);
    point = val(var_i);
end

% Get unique coordinates
x = cell(m_dim,1);
for i = 1 : m_dim
    x{i} = unique(obj.mesh(:,i));
end
% Get clustered mesh-points
y_data = get_var_dom(obj.local_idx, n_points, var_i);
% Works differently depending on some properties
if obj.is_mesh_variable && m_dim == 2
    % Plot original data
    try
        [xx, yy] = meshgrid(x{1}, x{2});
        F = scatteredInterpolant(obj.mesh(:,1), obj.mesh(:,2), y_data, 'natural');
        zz = F(xx, yy);
        figure();
        imshow(zz, [], 'XData', x{1}, 'YData', x{2}, 'Colormap', jet);
        axis xy; % Flip the figure
        axis on; 
    catch ME
        figure();
        scatter(obj.mesh(:,1), obj.mesh(:,2), 25, y_data, 'filled');
    end
    xlabel('x [m]'); ylabel('y [m]');
    title([arg_1, ' - Clustered mesh']);
elseif obj.is_mesh_variable && m_dim == 1
    figure();
    plot(obj.mesh, y_data, '--*', 'Color','k'); 
    grid on;
    xlabel('x [cm]','Interpreter','LaTex'); 
    ylabel('Cluster index [-]','Interpreter','LaTex');
    title(['\textbf{Clustered mesh: ', obj.variable_names{var_i},'}'],'Interpreter','LaTex');
    hold off;
elseif length(obj.local_idx) == size(obj.training_data, 2) && m_dim == 1
    % Plot original data
    figure();
    plot(x{1}, y_data, 'o');
    % Add stuff to the plot
    grid on;
    ylim([0, max(y_data)+1]);
    xlabel('x [m]'); ylabel('Cluster index [-]');
    title([num2str(point), ' - Clustered mesh ', var]);
else
    fprintf('Not possible to use this function if the mesh is not a variable.\n');
end

end

% Local functions
function y = get_var_dom(A, n, var_i)
y = A(1+n*(var_i-1):n*var_i, :);
end

function y = vectorize(X)

y = zeros(size(X,1) * size(X,2), 1);
for ii = 1 : size(X, 2)
    y(1+(ii-1)*size(X,1):(ii*size(X,1))) = X(:,ii);
end

end


