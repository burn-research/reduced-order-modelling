function [varargout] = get_variable_errors(obj, logical, var, model, varargin)

n_args = length(varargin);

% Provide this info
if n_args > 0
    is_mesh_var = varargin{1};
else
    is_mesh_var = true;
end

% Points
if obj.is_mesh_variable
    points_train = obj.training_points;
    points_pred = obj.prediction_points;
else
    points_train = unique( obj.training_points(:,size(obj.mesh,2)+1:end), 'rows' );
    points_pred = unique( obj.prediction_points(:,size(obj.mesh,2)+1:end), 'rows' );
end

% Get predicted/reconstructed variable
if logical
    Y = obj.get_variable(points_pred, var, model, is_mesh_var);
    Y_original = obj.get_variable(points_pred, var, 'data', is_mesh_var);
else
    Y = obj.get_variable(points_train, var, model, is_mesh_var);
    Y_original = obj.get_variable(points_train, var, 'data', is_mesh_var);
end

% Get error
varargout{1} = obj.getError(Y_original, Y, obj.mesh, var);

end


