function swapPoint(obj, varargin)
% Default values
points = [];

%% Input 
n_args = length(varargin);
if n_args == 0
    error('Not enough inputs.');
elseif n_args > 0
    points = varargin{1};
end

%% Main
imv = obj.is_mesh_variable;
[obj.training_points, obj.prediction_points, obj.training_data, obj.original_data]...
    = swapSample(points, obj.mesh, imv, obj.variable_names, ...
    obj.training_points, obj.prediction_points, obj.training_data, obj.original_data);

end


