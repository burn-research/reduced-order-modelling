function [y, varargout] = get_variable(obj, point, var, model, varargin)
%% Description
%{
    Get data of one specific physical variable, for a specifi point, from
    the correct data-set/model (original data, K-PCA model, K-CPCA model,
    etc).
%}

%% Function
if ~isempty(varargin)
    is_mesh_var = varargin{1};
else
    is_mesh_var = true;
end

% Find variable's name's position
var_i = 0;
for ii = 1 : length(obj.variable_names)
    if strcmp(var, obj.variable_names{ii})
        var_i = ii;
    end
end

% Mesh size
if iscell(obj.mesh)
    n = size(obj.mesh{var_i}, 1);
else
    n = size(obj.mesh, 1); 
end

% Understand which data-set must be used and get the correct column(s) I
if strcmp(model, 'pc')
    is_training = true;
    I = false(size(obj.pca_modes,2), 1);
    I(point) = true;
else
    dim = size(point,2); % Size of the input point
    % If mesh is variable, things are easier: look for 'point' with
    % ismember() and get I
    if obj.is_mesh_variable && (dim == size(obj.training_points, 2) || dim == size(obj.prediction_points, 2))
        [I, is_training] = get_I(point, obj.training_points, obj.prediction_points);
    elseif ~obj.is_mesh_variable && dim == (size(obj.training_points, 2) - size(obj.mesh,2))
        % If mesh is observations, we need to work differently depending on
        % the nature of the input 'point'. If the user has provided a
        % full point (comprised of mesh values): work as above. If the user
        % is looking only for the parameter values (with no mesh values),
        % some pre-processing is needed
        % Input point(s) with different dim: ismember() cannot work with
        % 'point' and the attributes 'training_points' and 'prediction_points'
        x_train = unique( obj.training_points(:,size(obj.mesh,2)+1:end), 'rows' );
        x_test = unique( obj.prediction_points(:,size(obj.mesh,2)+1:end), 'rows' );
        % Because of precision problems with ismember(), look for the
        % closest point
        [I, is_training] = get_I(point, x_train, x_test);
    elseif ~obj.is_mesh_variable && dim == size(obj.training_points, 2) 
        x = (point(:,size(obj.mesh,2)+1:end));
        x_train = unique( obj.training_points(:,size(obj.mesh,2)+1:end), 'rows' );
        x_test = unique( obj.prediction_points(:,size(obj.mesh,2)+1:end), 'rows' );
        [I, is_training] = get_I(x, x_train, x_test);
    else
        error('Should not have ended up here!');
    end
end

% Get the variable
switch model
    case 'centered'
        [A, B] = get_AB(is_mesh_var, obj.centered_scaled_data, obj.original_data, obj.mesh, obj.variable_names);
    case 'data'
        [A, B] = get_AB(is_mesh_var, obj.training_data, obj.original_data, obj.mesh, obj.variable_names);
    case 'pca'
        [A, B] = get_AB(is_mesh_var, obj.pca_recovered_data, obj.kpca_predictions, obj.mesh, obj.variable_names);
    case 'cpca'
        [A, B] = get_AB(is_mesh_var, obj.cpca_recovered_data, obj.kcpca_predictions, obj.mesh, obj.variable_names);
    case 'lpca'
        [A, B] = get_AB(is_mesh_var, obj.local_recovered_data_pca, obj.klpca_predictions, obj.mesh, obj.variable_names);
    case 'lcpca'
        [A, B] = get_AB(is_mesh_var, obj.local_recovered_data_cpca, obj.klcpca_predictions, obj.mesh, obj.variable_names);
    case 'direct'
        [A, B] = get_AB(is_mesh_var, [], obj.kriged_direct_data, obj.mesh, obj.variable_names);
    case 'pc'
        [A, B] = get_AB(is_mesh_var, obj.pca_modes, [], obj.mesh, obj.variable_names);
end

% In case A or B were empty, it means that the wanted model has not been
% run yet, and y will be an empty vector.

y = get_it(is_training, A, B, n, var_i, I);


end


function y = get_it(cond, A, B, n, var_i, I)

if cond
    if isempty(A)
        y = []; return;
    end
    iy1 = 1+n*(var_i-1);
    iy2 = n*var_i;
    y = A(iy1:iy2, I);
else
    if isempty(B)
        y = []; return;
    end
    iy1 = 1+n*(var_i-1);
    iy2 = n*var_i;
    y = B(iy1:iy2, I);
end

end


function [A, B] = get_AB(is_mesh_var, data_A, data_B, mesh, vars)

if is_mesh_var
    if isempty(data_A)
        A = [];
    else
        A = data_A;
    end
    if isempty(data_B)
        B = [];
    else
        B = data_B;
    end
else
    n_samples = size(data_A,2) / size(mesh,1);
    if isempty(data_A)
        A = [];
    else
        A = leaveStateVars(data_A, mesh, vars, n_samples);
    end
    n_samples = size(data_B,2) / size(mesh,1);
    if isempty(data_B)
        B = [];
    else
        B = leaveStateVars(data_B, mesh, vars, n_samples);
    end
end

end


function [I, is_training] = get_I(point, x_train, x_test)

% If the point is found, no problem
if ismember(point, x_train, 'rows')
    is_training = true;
    I = ismember(x_train, point, 'rows');
elseif ismember(point, x_test, 'rows')
    is_training = false;
    I = ismember(x_test, point, 'rows');
else
    % If the point is not found, look for the closest one
    I = [];
    [n_points, dim] = size(point);
    n_train = size(x_train, 1);
    n_test = size(x_test, 1);
    y = [x_train; x_test];
    % Do that for every point in 'point'
    for i = 1 : n_points
        s = (y - repmat(point(i,:), size(y,1), 1)).^2;
        dist = sqrt( sum(s, 2) ); % Distance of point(i) from all the others
        [~, idx] = min(dist); % Get the closest one
        % Understand if it was a training or prediction point
        if idx > n_train
            is_training = false;
            I = [I, idx - n_train]; % Get its position (offset)
        else
            is_training = true;
            I = [I, idx]; % Get its position (no offset)
        end
        % Small (wanted) bug: it is assumed that, in case 'point' is an array
        % of points, all of them are either training or prediction points.
        % If not, the info contained in 'is_training' is true only for the
        % last point(i) that's checked last
    end
    % We need logical indexing
    if is_training
        temp = false(n_train, 1);
        temp(I) = true;
    else
        temp = false(n_test, 1);
        temp(I) = true;
    end
    I = temp; % Here we go!
end

end





