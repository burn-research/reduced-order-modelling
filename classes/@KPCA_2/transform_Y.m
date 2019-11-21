function transform_Y(obj, varargin)

% Check input
if nargin > 0
    obj.xp_before_Q = varargin{1};
end
if nargin > 1
    obj.xp_kriged_before_Q = varargin{2};
end

% Retrieve important infos
n_samples_for_Y = size(obj.n_samples_for_Y, 1);
n_samples2_for_Y = size(obj.n_samples2_for_Y, 1);

% Get the Y matrix (Q, xq, vars, p)
obj.training_data = ...
    leaveStateVars(obj.training_data, obj.mesh, obj.variable_names, n_samples_for_Y);

% Get the Q_orig matrix
obj.original_data = ...
    leaveStateVars(obj.original_data, obj.mesh, obj.variable_names, n_samples2_for_Y);

%% Data
% Get the Q matrix
[obj.training_data, obj.training_points, xpnames] = ...
    leaveStateVars(obj.training_data, obj.mesh, obj.variable_names, obj.n_samples_for_Y);

% Get the Q_orig matrix
[obj.original_data, obj.prediction_points, xpnames_2] = ...
    leaveStateVars(obj.original_data, obj.mesh, obj.variable_names, obj.n_samples2_for_Y);

% Get Q for centered-scaled data
if ~isempty(obj.centered_scaled_data)
    obj.centered_scaled_data = ...
        leaveStateVars(obj.centered_scaled_data, obj.mesh, obj.variable_names, obj.n_samples_for_Y);
end

%% Direct Kriging
% Get Q for kriged_direct_data
if ~isempty(obj.kriged_direct_data)
    obj.kriged_direct_data = ...
        leaveStateVars(obj.kriged_direct_data, obj.mesh, obj.variable_names, obj.n_samples2_for_Y);
end

%% PCA
% Get Q for pca_recovered_data
if ~isempty(obj.pca_recovered_data)
    obj.pca_recovered_data = ...
        leaveStateVars(obj.pca_recovered_data, obj.mesh, obj.variable_names, obj.n_samples_for_Y);
end
% Get Q for kpca_predictions
if ~isempty(obj.kpca_predictions)
    obj.kpca_predictions = ...
        leaveStateVars(obj.kpca_predictions, obj.mesh, obj.variable_names, obj.n_samples2_for_Y);
end
% Local PCA
% Get Q for local_recovered_data_pca
if ~isempty(obj.local_recovered_data_pca)
    obj.local_recovered_data_pca = ...
        leaveStateVars(obj.local_recovered_data_pca, obj.mesh, obj.variable_names, obj.n_samples_for_Y);
    obj.local_recovered_centered_data = ...
        leaveStateVars(obj.local_recovered_centered_data, obj.mesh, obj.variable_names, obj.n_samples_for_Y);
end
% Get Q for klpca_predictions
if ~isempty(obj.klpca_predictions)
    obj.klpca_predictions = ...
        leaveStateVars(obj.klpca_predictions, obj.mesh, obj.variable_names, obj.n_samples2_for_Y);
end

%% CPCA
% Get Q for cpca_recovered_data
if ~isempty(obj.cpca_recovered_data)
    obj.cpca_recovered_data = ...
        leaveStateVars(obj.cpca_recovered_data, obj.mesh, obj.variable_names, obj.n_samples_for_Y);
end
% Get Q for kcpca_predictions
if ~isempty(obj.kcpca_predictions)
    obj.kcpca_predictions = ...
        leaveStateVars(obj.kcpca_predictions, obj.mesh, obj.variable_names, obj.n_samples2_for_Y);
end
% Local CPCA
% Get Q for local_recovered_data_cpca
if ~isempty(obj.local_recovered_data_cpca)
    obj.local_recovered_data_cpca = ...
        leaveStateVars(obj.local_recovered_data_cpca, obj.mesh, obj.variable_names, obj.n_samples_for_Y);
end
% Get Q for klcpca_predictions
if ~isempty(obj.klcpca_predictions)
    obj.klcpca_predictions = ...
        leaveStateVars(obj.klcpca_predictions, obj.mesh, obj.variable_names, obj.n_samples2_for_Y);
end

%% What about the errors?


%% Restore these quantities
obj.training_points = obj.xp_before_Q;
obj.prediction_points = obj.xp_kriged_before_Q;

%% Inform of this transformation
obj.is_mesh_variable = false;

end
