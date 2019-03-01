function getPcaErrors(obj, varargin)

[obj.pca_reconstruction_error_observations, ...
    obj.pca_reconstruction_error_variables] = obj.getError(...
    obj.training_data, obj.pca_recovered_data);

end


% OBSOLETE
% % This all if statement should become a funcion
% if obj.is_mesh_variable
%    % Get errors
%     obj.pca_reconstruction_error_observations = ...
%         obj.get_errors(obj.training_data, obj.pca_recovered_data, true); 
%     obj.pca_reconstruction_error_variables = ...
%         obj.get_errors(obj.training_data', obj.pca_recovered_data');
%     obj.pca_reconstruction_error_variables = ...
%         aveVar(obj.pca_reconstruction_error_variables, obj.variable_names, obj.mesh); 
% else
%     % Rebuild the Y matrices
%     n_samples = size(obj.training_data, 2) / size(obj.mesh,1);
%     Y_original = leaveStateVars(obj.training_data, obj.mesh, obj.variable_names, n_samples);
%     Y_predicted = leaveStateVars(obj.pca_recovered_data, obj.mesh, obj.variable_names, n_samples);
%     % Get errors
%     obj.pca_reconstruction_error_variables = ...
%         obj.get_errors(Y_original', Y_predicted');
%     obj.pca_reconstruction_error_observations = ...
%         obj.get_errors(Y_original, Y_predicted, true);
%     obj.pca_reconstruction_error_variables = ...
%         aveVar(obj.pca_reconstruction_error_variables, obj.variable_names, obj.mesh);
% end


