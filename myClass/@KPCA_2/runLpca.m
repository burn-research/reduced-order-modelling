function runLpca(obj, varargin)
%% Perform LPCA

% Warn that LPCA is being performed
fprintf('\nLocal PCA is being performed...\n');

% The property 'local_nz_idx_clust' will be re-evaluated now. This way, only 
% the property 'local_idx' is needed to run this subroutine.
if ~isempty(obj.local_idx)
    k = length(unique(obj.local_idx));
else
    k = 1;
end
if k > 1
    obj.local_nz_idx_clust = cell(k, 1);
    for i = 1 : k 
        obj.local_nz_idx_clust{i} = find(obj.local_idx == i);
    end
else
    obj.local_nz_idx_clust = [];
end

% Build the Local PCA objects
if ~isempty(obj.local_nz_idx_clust)
    % OPTIONS is a structure, needed to pass user-provided mean_column and
    % scaling_factors. We do not want to center-scale data in the local
    % clusters, but then need external mean and scaling factors if we wish
    % to rebuild data or apply CPCA in a cluster. CLUSTERS is a cell-array
    % of training data.
    
    % Need to find a way to understand the number of clusters that is going
    % to be used now for LPCA
    current_numberOfClusters = length(obj.local_nz_idx_clust);
    obj.number_of_clusters = current_numberOfClusters;
    
    % Initialize this variables accordingly
    options = cell(current_numberOfClusters, 1);
    clusters = cell(current_numberOfClusters, 1);
    
    % Set the options for each cluster 
    for i = 1 : current_numberOfClusters
        % Get the data (set of rows or columns)
        clusters{i} = obj.get_centered_cluster(i); % obj.get_cluster(i);
        % The mean vector and scaling factors need different operations
        % w.r.t. the clustering dimension. They are passed mostly because
        % in case of Local CPCA, it takes to know info about the real
        % variables.
        % TO DO: do not pass them now, only in runLcpca()
        if ~obj.clustering_dimension
            % Pass a subset of rows of the mean vector and scaling factors
            options{i}.mean_column = obj.mean_column(obj.local_nz_idx_clust{i});
            options{i}.scaling_factors = obj.scaling_factors(obj.local_nz_idx_clust{i});
        else
            % Pass all of it
            options{i}.mean_column = obj.mean_column;
            options{i}.scaling_factors = obj.scaling_factors;
        end
        % Set the created object as Local PCA and pass info about the mesh
        options{i}.is_local = true;
        options{i}.is_mesh_variable = obj.is_mesh_variable;
    end
    
    % Create the object in parallel for speed
    obj.local_pca = parallel_lpca(current_numberOfClusters, ...
        obj.clustering_dimension, ...
        options, ...
        clusters, ...
        obj.pca_approximation_order, ...
        obj.trend_function, ...
        obj.correlation_function, ...
        obj.correlation_function_optimization, ...
        obj.centered_scaled_data, ...
        obj.training_points, ...
        obj.variable_names, ...
        size(obj.mesh, 1), ...
        obj.is_mesh_variable, ...
        obj.scaling_criterion, ...
        obj.hyp, ...
        obj.meanfunc, ...
        obj.covfunc, ...
        obj.likfunc, ...
        obj.inffunc, ...
        obj.hyp_guess);
    
    % In case the test data is there and rows have been clustered, pass the
    % corresponding rows of the test data (why do this at all?)
%     if ~isempty(obj.original_data) && ~obj.clustering_dimension
%         obj.local_pca{i}.original_data = obj.original_data(obj.local_nz_idx_clust{i},:);
%     end
    
else
    obj.local_pca = {};
    warning('No clusters were found, try to runCLustering() first.');
    return;
end

fprintf('\nLocal PCA has been performed.\n');

%% Update some properties
obj.update_lpca_prop();

%% Recover data
obj.recoverLpca();

%% Estimate errors
try
    obj.getLpcaErrors();
catch ME
    fprintf('[getLpcaErrors:] Could not getPcaErrors()');
    disp(ME)
    for ii = 1 : length(ME.stack)
        disp(ME.stack(ii));
    end
end

end




function lpca = parallel_lpca(num_clusts, clustering_dimension, ...
    options, ...
    clusters, k, trend_function, correlation_function, ...
    correlation_function_optimization, Ycs, xp, vars, n_points, is_mesh_var, sc, ...
    hyp, meanfunc, covfunc, likfunc, inffunc, hyp_guess)

lpca = cell(num_clusts, 1);

parfor i = 1 : num_clusts
    % LPCA
    lpca{i} = KPCA_2(clusters{i}, options{i}); % Create the LocalPCA objects
    lpca{i}.pca_approximation_order = k;
    lpca{i}.trend_function = trend_function; % Set the same trend function
    lpca{i}.correlation_function = correlation_function; % Set the same correlation function
    lpca{i}.correlation_function_optimization = correlation_function_optimization;
    lpca{i}.center_scale = false;
    lpca{i}.variable_names = vars;
    lpca{i}.mesh = n_points;
    lpca{i}.is_mesh_variable = is_mesh_var;
    lpca{i}.is_local = true;
    lpca{i}.scaling_criterion = sc;
    lpca{i}.hyp = hyp;
    lpca{i}.meanfunc = meanfunc;
    lpca{i}.covfunc = covfunc;
    lpca{i}.likfunc = likfunc;
    lpca{i}.inffunc = inffunc;
    lpca{i}.hyp_guess = hyp_guess;
    
    % Run PCA
    lpca{i}.runPca();
    
    % Useful for the Q-matrix
    if clustering_dimension
        lpca{i}.morevalues = lpca{i}.pca_modes' * Ycs;
        lpca{i}.moresamples = xp;
    end       
end

end



