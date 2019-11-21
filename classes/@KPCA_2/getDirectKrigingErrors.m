function getDirectKrigingErrors(obj, varargin)

% If there is no test data, return
if isempty(obj.original_data)
    return
end

if obj.is_mesh_variable
    % Get errors
    obj.direct_prediction_error_observations = ...
        obj.get_errors(obj.original_data, obj.kriged_direct_data, true);
    obj.direct_prediction_error_variables = ...
        obj.get_errors(obj.original_data', obj.kriged_direct_data');
    obj.direct_prediction_error_variables = ...
        aveVar(obj.direct_prediction_error_variables, obj.variable_names, obj.mesh);
else
    % Rebuild the Y matrices
    n_samples = size(obj.kpca_predictions,2) / size(obj.mesh,1);
    Y_original = leaveStateVars(obj.original_data, obj.mesh, obj.variable_names, n_samples);
    Y_predicted = leaveStateVars(obj.kriged_direct_data, obj.mesh, obj.variable_names, n_samples);
    % Get errors
    obj.direct_prediction_error_observations = obj.get_errors(Y_original, Y_predicted, true);
    obj.direct_prediction_error_variables = obj.get_errors(Y_original', Y_predicted');
    obj.direct_prediction_error_variables = ...
        aveVar(obj.direct_prediction_error_variables, obj.variable_names, obj.mesh);
end

end

