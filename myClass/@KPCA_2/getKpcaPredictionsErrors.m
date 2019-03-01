function getKpcaPredictionsErrors(obj, varargin)

% If there is no test data, return
if isempty(obj.original_data)
    return
end

[obj.kpca_prediction_error_observations, ...
    obj.kpca_prediction_error_variables] = obj.getError(...
    obj.original_data, obj.kpca_predictions);

end


% OBSOLETE
% % This all if statement should become a funcion
% if obj.is_mesh_variable
%    % Get errors
%     obj.pca_reconstruction_error_observations = ...
%         obj.get_errors(obj.original_data, obj.kpca_predictions, true); 
%     obj.pca_reconstruction_error_variables = ...
%         obj.get_errors(obj.original_data', obj.kpca_predictions');
%     obj.pca_reconstruction_error_variables = ...
%         aveVar(obj.pca_reconstruction_error_variables, obj.variable_names, obj.mesh); 
% else
%     % Rebuild the Y matrices
%     n_samples = size(obj.original_data, 2) / size(obj.mesh,1);
%     Y_original = leaveStateVars(obj.original_data, obj.mesh, obj.variable_names, n_samples);
%     Y_predicted = leaveStateVars(obj.kpca_predictions, obj.mesh, obj.variable_names, n_samples);
%     % Get errors
%     obj.pca_reconstruction_error_variables = ...
%         obj.get_errors(Y_original', Y_predicted');
%     obj.pca_reconstruction_error_observations = ...
%         obj.get_errors(Y_original, Y_predicted, true);
%     obj.pca_reconstruction_error_variables = ...
%         aveVar(obj.pca_reconstruction_error_variables, obj.variable_names, obj.mesh);
% end


