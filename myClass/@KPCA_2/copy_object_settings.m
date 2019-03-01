function [varargout] = copy_object_settings(obj, varargin)
%% Input
n_args = length(varargin);
if n_args > 0 && ~isempty(varargin{1})
    kpca = varargin{1};
else
    kpca = KPCA_2([]);
end
if n_args > 1
    copy_data_info = varargin{2};
else
    copy_data_info = false;
end

%% Copying properties
% Copy mesh, training and prediction points
if copy_data_info
    kpca.mesh = obj.mesh;
    kpca.training_points = obj.training_points;
    kpca.prediction_points = obj.prediction_points;
    kpca.original_data = obj.original_data;
end

% Input/Output space
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

%% Output
if nargout > 0
    varargout{1} = kpca;
end

end


