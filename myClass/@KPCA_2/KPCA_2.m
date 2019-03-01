classdef KPCA_2 < handle
%% Description


%% Properties
properties (Constant)
    % SM models
    sm_models = {'pca', 'cpca', 'lpca', 'lcpca', 'direct'};
    
    % Regression methods
    reg_m = {{'Kriging', 'kriging', 'krig', 'k'};...
             {'Linear', 'linear', 'MARS', 'mars'};...
             {'gpml', 'GPML', 'Gpml', 'g'}};
    
    % Kriging Toolbox
    correlation_functions = {@correxp; @corrgauss; @corrmatern32; ...
        @corrmatern52; @corrmatern12; @corrlin; @corrspline};
    trend_functions = {'regpoly0', 'regpoly1', 'regpoly2', 'regpoly3', 'regpoly4'};
    
    % GPML Toolbox
    covariance_functions = {};
    mean_functions = {@meanZero, @meanOne, @meanConst, @meanLinear,...
        @meanPoly, @meanDiscrete, @meanGP, @meanGPexact, @meanNN, @meanWSPC};
end

properties (Access = public)
    % Numerical
    scaled_errors = false;
    r_tol = 1e-9; % Used for evaluation of errors
    is_parallel = true;
    
    % Input/Output space
    training_points = []; % Input values for the scores
    prediction_points = []; % Input values for the scores to be Kriged 
    mesh = 1; % Mesh
    original_data = []; % Original data: true solutions (no repeated data)
    variable_names = {}; % Names of the physical variables
    parameter_names = {}; % Names of the physical parameters
    is_mesh_variable = true; % Spatial mesh status
    
    % PCA
    pca_manual = false;
    center_scale = true;
    scaling_criterion = 1;
    pca_approximation_order = 1; % PCA approximation order
%     log_transform = false; % use a log-transformation to keep data positive
    
    % Regression/Interpolation
    reg_method = 'Kriging';
    is_kriging = true;
    is_gpml = false;
    is_linear = false;
    stacking = false;
    includeMoreSamples = false; % Include more samples to improve interpolation/regression
    targetBYtarget = false;
    cs_samples = true;
    selected_mesh_points = [];
    cluster_prediction = 0; % 0) Closest; 1) 2) 3)
    
    % Kriging
    hyp_guess = []; % If empty, it will be generated automatically [check generateHyperparameters0()]
    is_interpolation = true;
    trend_function = 'regpoly0'; % Kriging: trend function
    correlation_function = @correxp; % Kriging: correlation function
    correlation_function_optimization = false; % Correlation Function Optimization    
    
    % MARS
    mars_maxOrd = 3;
    
    % GPML
    meanfunc = {@meanSum, {@meanLinear, @meanConst}}; % Mean function
    covfunc = @covSEard; % Covariance function
    likfunc = @likGauss; % Likelihood
    inffunc = @infExact; 
    hyp = struct('mean', [], 'cov', [], 'lik', -1);
    h_optimize = false;
    
    % Constrained PCA
    my_constraint = @KPCA_2.allPositiveCon; % Constraint function
    cpca_guess_corrector = 1; % Description: (guess = corrector * obj.pca_scores)
    cpca_all = false; % Evaluate the CPCA scores for all values of the app. order
    
    % Clustering
    clustering_dimension = false; % Clustering on variables or samples
    number_of_clusters = 1; % Number of clusters
    accelerated_clustering = false;
    
    % Local PCA
    centroid = [];
    idx_0 = []; % Initial indeces for Local PCA
    lpca_clustering_options = false; % Stop at kmeans(): false, true
    dist_on_cs_points = false;
end

properties (SetAccess = public)
    % Numerical
    training_data; % Training data
    centered_scaled_data; % Centered and scaled data
    mean_column; % Mean column
    scaling_factors; % Scaling factors
    correlation_matrix; % Correlation matrix
    
    % PCA
    pca_eigenvalues; % Eigenvalues
    pca_modes; % PCA modes
    pca_scores; % PCA scores
    pca_scores_kriging_targets;
    pca_energy; % Captured variance
    pca_recovered_data; % Recovered data from PCA
    pca_reconstruction_error_observations;  %
    pca_reconstruction_error_variables; %

    % Constrained PCA
    cpca_scores; % CPCA scores
    cpca_scores_stored = {};
    cpca_recovered_data; % Recovered data from CPCA
    cpca_reconstruction_error_observations; % 
    cpca_reconstruction_error_variables; %
    
    % Local PCA
    local_idx = []; % Cluster assignments
    local_idx_stored = {}; % Set of cluster assignments
    number_of_clusters_stored = [];
    local_nz_idx_clust; % Cell array of clustered rows
    local_pca; % Cell array of KPCA objects
    local_recovered_centered_data;
    local_recovered_data_pca; %
    local_recovered_data_cpca; %
    local_pca_reconstruction_error_observations; %
    local_pca_reconstruction_error_variables; %
    local_cpca_reconstruction_error_observations; %
    local_cpca_reconstruction_error_variables; %
    
    % Kriging
    kriging_weigths = []; % Weights for different predictions/regression functions
    pcaKrigingModel;
    pcaKrigingMSE;
    cpcaKrigingModel;
    cpcaKrigingMSE;
    directKrigingModel;
    directKrigingMSE;
    kriged_pca_scores; % Kriged PCA scores
    kriged_cpca_scores; % Kriged CPCA scores
    mse; % Variance of the Kriging interpolation
    kriged_direct_data; % Kriging applied directly on the original data
    direct_prediction_error_observations;
    direct_prediction_error_variables;
    kpca_predictions; %
    kpca_prediction_error_observations; %
    kpca_prediction_error_variables;
    kcpca_predictions; %
    kcpca_prediction_error_observations; %
    kcpca_prediction_error_variables; %
    klpca_predictions; %
    klpca_prediction_error_observations; % 
    klpca_prediction_error_variables; %
    klcpca_predictions; %
    klcpca_prediction_error_observations; %
    klcpca_prediction_error_variables; %
    
    % Cross-Validation
    pcaCV;
    cpcaCV;
    directCV;
    lpcaCV;
    lcpcaCV;
    
    % Rebuilt and sorted data
    pca_scores_sorted; % Sorted PCA scores and Kriged PCA scores
    sorted_points; % Sorted training and prediction points
    sorted_kpca_data; % Sorted PCA recovered and Kriged data
    sorted_kriged_direct_data; % Sorted data, Y_k_direct
    sorted_klpca_data; % Full reconstruction of the data from the K+LPCA
    sorted_kcpca_data; 
end


%% Additional properties
properties (SetAccess = private, Hidden)
    is_local = false; % Is this a local cluster or not?
    
    % Needed for the Kriging to work better
    moresamples = [];           % Location of the more samples
    morevalues = [];            % Values of the more samples
    
    % Some saved bytes
    xp_before_Q;
    xp_kriged_before_Q;
    xp_before_Y;
    xp_kriged_before_Y;
end


%% Constructor
methods
    function obj = KPCA_2(data, varargin)
        % Load data
        obj.training_data = data;
        % Check for further inputs
        if ~isempty(varargin)
            options = varargin{1}; % Options structure
            % Column of mean scalars
            if isfield(options, 'mean_column')
                obj.mean_column = options.mean_column;
                obj.center_scale = false;
            end
            % Vector of scaling factors
            if isfield(options, 'scaling_factors')
                obj.scaling_factors = options.scaling_factors;
                obj.center_scale = false;
            end
            if isfield(options, 'is_mesh_variable')
                obj.is_mesh_variable = options.is_mesh_variable;
            end
            % Logical flag
            if isfield(options, 'is_local')
                obj.is_local = options.is_local;
                obj.center_scale = false;
            end
        end
        obj.setMaxApproxOrder(); % Default is 1, thus get the maximum
    end
end


%% Set methods
methods
    % GPML
    % Meanfunc
    function set.meanfunc(obj, val)
        if isnumeric(val)
            if val == 0
                obj.meanfunc = @meanConst;
            elseif val == 1
                obj.meanfunc = {@meanSum, {@meanLinear, @meanConst}};
            elseif val == 2
                obj.meanfunc = {@meanSum, {@meanPoly, @meanConst}};
            end
        else
            obj.meanfunc = val;
        end
    end
    
    % Scaling criterion
    function set.scaling_criterion(obj, val)
        % Different choices are available:
        % -1) Use the mean_column as mean vector, whatever is there already
        % 0) No scaling 
        % 1) Auto-scaling (STD), each variable is normalized by its standard
        % deviation 
        % 2) RANGE each variable is normalized by its range 
        % 3) PARETO, each variable is scaled by the square root of its standard
        % deviation  
        % 4) VAST, each variable is scaled by the standard deviation and
        % coefficient of variation 
        % 5) LEVEL, each variable is normalized by the mean of the data
        % 6) MAX, each variable is scaled by its maximum value
        if val ~= round(val) || val < -1 || val > 9
            error('Unknown scaling criterion');
        end
        obj.scaling_criterion = val;
    end
    
    % I/O: training points
    function set.training_points(obj, val)
        if isempty(val)
            obj.training_points = val;
            return
        end
%         [m, ~] = size(val);
%         [~, n] = size(obj.training_data);
%         if m ~= n
%             error('Number of rows of this input should equal the number of columns of the training data matrix.');
%         end
        obj.training_points = val;
    end
    
    % I/O: prediction points
    function set.prediction_points(obj, val)
        if isempty(val)
            obj.prediction_points = val;
            return
        end
%         [~, m] = size(val);
%         [~, n] = size(obj.training_points);
%         if m ~= n
%             error('The training and prediction points have different dimensions.');
%         end
        obj.prediction_points = val;
    end
    
    % I/O: mesh
    function set.is_mesh_variable(obj, val)
        if ~islogical(val) && val ~= 0 && val ~= 1
            error('This property must be either TRUE or FALSE.');
        end
        obj.is_mesh_variable = val;
    end
    
    % I/O: original data
    function set.original_data(obj, val)
        if isempty(val) || obj.is_local
            obj.original_data = val;
            return
        end
%         [r, c] = size(val);
%         [r_prediction_points, ~] = size(obj.prediction_points);
%         [r_training_data, ~] = size(obj.training_data);
%         if r ~= r_training_data || c ~= r_prediction_points
%             error('Input has unacceptable size.');
%         end
        obj.original_data = val;
    end
    
    % PCA: approximation order
    function set.pca_approximation_order(obj, f)
        if size(obj.training_data,1) > 0 && size(obj.training_data,2) > 0
            max_app = size(obj.training_data, 2);
        else
            max_app = 1;
        end
        if size(obj.pca_modes,1) > 0 && size(obj.pca_modes,2) > 0
            max_app = size(obj.pca_modes, 2);
        end
        if ~isnumeric(f)
            warning('Must provide a numeric input. Set to 1.');
            obj.pca_approximation_order = max([max_app, 1]);
        end
        if f >= 1
            obj.pca_approximation_order = f;
        elseif f <= 0
            warning('Pca approximation order must be > 0. Set to 1.');
            obj.pca_approximation_order = max_app;
        else
            obj.pca_approximation_order = obj.getApproximationOrder(f);
        end
    end 
    
    % PCA: center_scale
    function set.center_scale(obj, f)
        if ~islogical(f)
            warning('The property CENTER_SCALE must be a LOGICAL. Left unchanged.');
            return;
        end
        obj.center_scale = f;
    end 
    
    % Kriging: correlation function
    function set.correlation_function(obj, f)
        if ~isa(f,'char') && ~isa(f,'function_handle')
            f = obj.correlation_functions{f};
        end
        this = false(length(obj.correlation_functions), 1);
        for i = 1 : length(obj.correlation_functions)
            this(i) = isequal(f, obj.correlation_functions{i});
        end
        if ~any(this)
            fprintf('No possiblen correlation function provided.\n\n');
            return;
        end
        obj.correlation_function = f;
    end
    
    % Kriging: trend function
    function set.trend_function(obj, f)
        if isa(f, 'double')
            obj.trend_function = obj.trend_functions{f+1};
        else
            obj.trend_function = 'regpoly0';
        end
    end
    
    % Kriging: correlation function optimization
    function set.correlation_function_optimization(obj, val)
        if ~islogical(val) && val ~= 0 && val ~= 1
            error('This property must be either TRUE or FALSE.');
        end
        obj.correlation_function_optimization = val;
    end
    
    % Local PCA: clustering dimension
    function set.clustering_dimension(obj, val)
        if ~islogical(val) && val ~= 0 && val ~= 1
            error('This property must be either TRUE or FALSE.');
        end
        obj.clustering_dimension = val;
    end
    
    % Mesh & variable_names
    function set.mesh(obj, val)
        if obj.is_local
            obj.mesh = val;
            return
        end
        if isempty(obj.variable_names)
            obj.mesh = val;
            return
        elseif iscell(val)
            obj.mesh = val;
        elseif obj.is_mesh_variable && (size(obj.training_data,1) / size(val,1)) ~= length(obj.variable_names)
%             error('[@KPCA_2:set.mesh] There seems to be something wrong.');
        end
        obj.mesh = val;
    end
%     function set.variable_names(obj, val)
%         if isempty(obj.mesh)
%             obj.variable_names = val;
%             return
%         elseif (size(obj.training_data,1) / size(obj.mesh,1)) ~= length(val)
%             error('[@KPCA_2:set.variable_names] There seems to be something wrong.');
%         end
%         obj.variable_names = val;
%     end
end


%% Get methods
methods

end


%% Methods
methods
    % Data transforms
    [varargout] = transformData(obj, fun, varargin);
    
    % General routines
    runReconstruction(obj, varargin);
    runAll(obj, varargin);
    runAllConstrained(obj, varargin);
    setIdx(obj, idx);
    setClusterAssignments(obj, nc_int);
    setMaxApproxOrder(obj, varargin);
    setScores(obj, varargin);
    setCpcaScores(obj, varargin);
    [varargout] = applyConstraint(obj, varargin);
    eliminateVariable(obj, varargin);
    swapPoint(obj, varargin);
    setRegressionMethod(obj, method, varargin);
    getPredictionErrors(obj, varargin);
    [varargout] = copy_object_settings(obj, varargin);
    
    % PCA routines
    runPca(obj, varargin);
    runCpca(obj, varargin);
    runClustering(obj, varargin);
    runLpca(obj, varargin);
    runLcpca(obj, varargin);
    cleanPca(obj, varargin);
    cleanCpca(obj, varargin);
    cleanLpca(obj, varargin);
    cleanLcpca(obj, varargin);
    
    % Kriging routines
    runPcaKriging(obj, varargin);
    runCpcaKriging(obj, varargin);
    runLpcaKriging(obj, varargin);
    runLcpcaKriging(obj, varargin);
    runDirectKriging(obj, varargin);
    get_Kriging_targets(obj, varargin);
    CVPE(obj, model, varargin);
    
    % Kriging predict only
    [varargout] = predict_PCA(obj, varargin);
    [varargout] = predict_LPCA(obj, varargin);
    [varargout] = predict_CPCA(obj, varargin);
    [varargout] = predict_LCPCA(obj, varargin);
%     trainPcaKriging(obj, varargin);
%     trainCpcaKriging(obj, varargin);
%     trainLpcaKriging(obj, varargin);
%     trainLcpcaKriging(obj, varargin);
%     trainDirectKriging(obj, varargin);
    
    % HDF5 routines
    write(obj, data, varargin);

    % Plot routines
    plot_captured_variance(obj, varargin);
    plot_eigenspectrum(obj, varargin);
    plot_compare_variable(obj, var, point, varargin);
    plot_sampled_input_space(obj, varargin);
    plot_parity(obj, char_var, boolean_var, varargin);
    plot_Kriging_PCA_scores(obj, variable_name, point, varargin);
    [varargout] = plot_contour(obj, varargin);
    [varargout] = plot_predictionErrors(obj, varargin);
    [varargout] = plot_clustered_vardomain(obj, varargin);
    [varargout] = plot_reconstructionErrors(obj, varargin);
    plot_cluster_populations(obj, varargin);
    [varargout] = plot_clustered_variable_space(obj, char_var, varargin);
    spy_centered_data(obj, varargin);
    [varargout] = plot_pca_mode(obj, mode_idx, var, varargin);
    
    % Useful routine
    [varargout] = get_variable_errors(obj, var, data_string, varargin);
end

methods (Hidden)
    [varargout] = getError(obj, varargin);
    k = train(obj, i, varargin);
    setPrivateProperty(obj, property_name, value);
    runCenterScale(obj, varargin);
    this = uncenterUnscale(obj, varargin);
    [decoded_data, varargout] = pca_decode(obj, modes, scores, scaling_factors, mean_column, varargin);
    val = getApproximationOrder(obj, f);
    update_gamma(obj, varargin);
    [varargout] = update_Kriging(obj, first_input, second_input, varargin);
    getKrigingPoints(obj);
    [idx, uncentered_nz_X_k, varargout] = localPCA_hierarchical(X, n_eigs, varargin);
    localPCA_accelerated(obj, varargin);
    choosePartitioningCriteria(obj, X, n_eigs, varargin);
    out = create_lpca_object(X, options);
    runPca_hiddenFunction(obj, varargin);
    recoverPca(obj, varargin);
    recoverCpca(obj, varargin);
    recoverLpca(obj, varargin);
    [cluster, varargout] = get_cluster(obj, cluster_index);
    [cluster, varargout] = get_centered_cluster(obj, cluster_index);
    getKpcaPredictions(obj, varargin);
    getKpcaPredictionsErrors(obj, varargin);
    getKlpcaPredictions(obj, varargin);
    getKcpcaPredictions(obj, varargin);
    getKlcpcaPredictions(obj, varargin);
    getKcpcaPredictionsErrors(obj, varargin);
    getDirectKrigingErrors(obj, varargin)
    y = get_errors(obj, x_original, x_predicted, varargin);
    y = get_error(obj, x_original, x_predicted, varargin);
    getPcaErrors(obj, varargin);
    getCpcaErrors(obj, varargin);
    getLpcaErrors(obj, varargin);
    getLcpcaErrors(obj, varargin);
    plot_reconstruction_error(obj, y1, y2, title_string);
    plot_pca_reconstruction_error(obj, varargin);
    plot_cpca_reconstruction_error(obj, varargin);
    plot_lpca_reconstruction_error(obj, varargin);
    plot_lcpca_reconstruction_error(obj, varargin);
    [y, varargout] = get_variable(obj, point, var, data_string, varargin);
    centroid = getCentroid(obj, varargin);
    centered_scaled_data = subtractCentroid(obj, varargin);
    [varargout] = plot_inputSpaceErrors(obj, varargin);
    setAllPositive(obj, varargin);
    [varargout] = getClusterPopulations(obj, varargin);
end

methods (Static)
    [c, ceq] = allPositiveCon(gamma, obj, varargin);
    [predictions, mse, k] = linearRegression(samples, values, prediction_points, varargin);
    [krigingOutput, mse, k] = stackKriging(samples, values, p_points, t_fun, c_funs, guesses, varargin);
    [krigingOutput, mse, k] = ensembleKriging(samples, values, prediction_points, trendFuns, corrFuns, guess, varargin);
end


end

