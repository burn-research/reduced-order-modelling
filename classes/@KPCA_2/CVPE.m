function CVPE(obj, model, varargin)
%% Description
% This routine performs Cross-Validation. A surrogate model is built from 
% N-1 training observations. The N-th observation is used as test data. The
% process needs to be repeated N times.
%

%% Pre-processing
imv = obj.is_mesh_variable;
n_vars = size(obj.training_data,1);
[m_points, dim_mesh] = size(obj.mesh); % Number of mesh points
% Get general form for the training points (get rid of the mesh)
if imv
    training_points = obj.training_points;
else
    training_points = unique(obj.training_points(:,dim_mesh+1:end), 'rows');
end
% Get number of total training observations
N = size(training_points,1);
err_collect = zeros(n_vars,N);

%% Main (PCA calculations are done again)
for ii = 1 : N
    kpca = KPCA_2([]); % Create KPCA_2 object with no data
    % Set KPCA propertites
    kpca = copy_object(kpca, obj);
    % Leave one out
    [kpca.training_points, kpca.prediction_points, kpca.training_data, kpca.original_data]...
        = swapSample(training_points(ii,:), obj.mesh, imv, obj.variable_names, ...
        obj.training_points, [], obj.training_data, []);
    % NOTES: swapSamples() sticks to the format set by imv
    % Get the difference between the test data and the prediction. In case
    % of imv == true, a matrix (not a vector) is return, hence the mean().
    % A matrix is returned becasue one point in that case means one
    % parameter point but the whole mesh.
    fprintf('\nLeaving out %d (out of %d).', ii, N);
    switch model
        case 'pca'
            evalc('kpca.runPca(); kpca.runPcaKriging()');
            err_collect(:,ii) = mean(kpca.original_data - kpca.kpca_predictions, 2);
        case 'cpca'
            % [Missing]
        case 'lpca'
            evalc('kpca.runLpca(); kpca.runLpcaKriging()');
            err_collect(:,ii) = mean(kpca.original_data -kpca.klpca_predictions, 2);
        case 'lcpca'
            % [Missing]
        case 'direct'
            evalc('kpca.runDirectKriging();');
            err_collect(:,ii) = mean(kpca.original_data - kpca.kriged_direct_data, 2);
    end
    kpca = [];
end

%% Main (no PCA calculations, just leave one set of scores out)


%% Post-processing
switch model
    case 'pca'
        obj.pcaCV = NRMSE(err_collect, obj.scaling_factors);
    case 'cpca'
        % [Missing]
    case 'lpca'
        obj.lpcaCV = NRMSE(err_collect, obj.scaling_factors);
    case 'lcpca'
        % [Missing]
    case 'direct'
        obj.directCV = NRMSE(err_collect, obj.scaling_factors);
end
fprintf('\n');

end

function y_err = NRMSE(delta, sigma)
% error = sqrt( 1/s * sum_i^s[ (rom - full)^2 ] ) / (max_full - min_full)
a_tol = 1e-16;
[s, ~] = size(delta);
% Get errors
y_err = zeros(s,1);
for i = 1 : s
    y_err(i) = sqrt( (1/s) * sum(delta(i,:).^2) ) / (sigma(i) + a_tol);
end
y_err(y_err < a_tol) = 0;
end

function kpca = copy_object(kpca, obj)
% Input/Output space
kpca.mesh = obj.mesh;
kpca.is_mesh_variable = obj.is_mesh_variable;
kpca.variable_names = obj.variable_names;

% Regression method
kpca.setRegressionMethod(obj.reg_method);

% PCA
kpca.pca_manual = obj.pca_manual;
kpca.center_scale = obj.center_scale;
kpca.scaling_criterion = obj.scaling_criterion;
kpca.pca_approximation_order = obj.pca_approximation_order;

% Regression
kpca.cs_samples = obj.cs_samples;
kpca.selected_mesh_points = obj.selected_mesh_points;
kpca.cluster_prediction = obj.cluster_prediction; 

% Kriging
kpca.hyp_guess = kpca.hyp_guess;
kpca.is_interpolation = obj.is_interpolation;
kpca.trend_function = obj.trend_function; % Kriging: trend function
kpca.correlation_function = obj.correlation_function; % Kriging: correlation function
kpca.correlation_function_optimization = obj.correlation_function_optimization; % Correlation Function Optimization    

% MARS
kpca.mars_maxOrd = obj.mars_maxOrd;

% GPML
kpca.meanfunc = obj.meanfunc;
kpca.covfunc = obj.covfunc;
kpca.likfunc = obj.likfunc;
kpca.inffunc = obj.inffunc;
kpca.hyp = obj.hyp;

% Constrained PCA
kpca.my_constraint = @KPCA_2.allPositiveCon; % Constraint function
kpca.cpca_guess_corrector = 1; % Description: (guess = corrector * obj.pca_scores)
kpca.cpca_all = false;

% Clustering
kpca.clustering_dimension = obj.clustering_dimension; % Clustering on variables or samples
kpca.number_of_clusters = obj.number_of_clusters; % Number of clusters
kpca.accelerated_clustering = obj.accelerated_clustering;

% Local PCA
kpca.centroid = obj.centroid;
kpca.idx_0 = obj.idx_0; % Initial indeces for Local PCA
kpca.lpca_clustering_options = obj.lpca_clustering_options;
kpca.dist_on_cs_points = obj.dist_on_cs_points;

end







