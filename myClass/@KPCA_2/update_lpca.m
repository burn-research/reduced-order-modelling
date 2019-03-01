function obj = update_lpca(obj, varargin)

lpca_metric = obj.lpca_metric;  % Metric for from-cluster distances
pcascores_space = false;        % Perform partition in the PCA scores' space


if nargin > 1
    lpca_metric = varargin{1};
    if nargin > 2
        pcascores_space = varargin{2};
    end
end


% Check the number of clusters
if obj.num_clusts == 1
    obj.idx = ones(size(obj.Y, obj.clustering_dim), 1);
    obj.clusters = {};
    obj.lpca = {};
    obj.nz_idx_clust = {};
    return; % Do not run the clustering optimization procedure
end


% What follows runs only if NUM_CLUSTS > 1
%% Clustering optimization procedure

if obj.clustering_dim == 1 && obj.lpca_clustering_options == 0
    % Partition variables
    [obj.idx, obj.clusters, this] = KPCA.localPCA(obj.Ycs, obj.k, lpca_metric, pcascores_space, obj.num_clusts, obj.idx_0);

elseif obj.clustering_dim == 2 && obj.lpca_clustering_options == 0
    % Partion observations
    [obj.idx, obj.clusters, this] = KPCA.localPCA(obj.Ycs', obj.k, lpca_metric, pcascores_space, obj.num_clusts, obj.idx_0);

elseif obj.lpca_clustering_options == 1 && obj.x_status_double == 1
    % Partition physical variables
    Q = getStateVars(obj.Ycs, obj.xq, obj.vars, obj.cols);
    [idx, clusters, this] = KPCA.localPCA(Q, obj.k, lpca_metric, pcascores_space, obj.num_clusts, obj.idx_0);
    for i = 1 : length(clusters)
        clusters{i} = leaveStateVars(clusters{i}, obj.xq, obj.vars, obj.cols);
        temp = [];
        for j = 1 : length(this.nz_idx_clust{i})
            temp = [temp, repmat(this.nz_idx_clust{i}(j), length(obj.xq), 1)];
        end
        this.nz_idx_clust{i} = temp;
    end
    obj.clusters = clusters;
    temp = [];
    for i = 1 : length(idx)
        temp = [temp, repmat(idx(i), length(obj.xq), 1)];
    end
    obj.idx = temp;

else 
    error('Property CLUSTERING_DIM was found being neither 1 nor 2, or LPCA_CLUSTERING_OPTIONS was 1 but X_STATUS is not VARIABLE.'); 
end

obj.nz_idx_clust = this.nz_idx_clust;


% Check the number of clusters returned is the same
if length(obj.clusters) ~= obj.num_clusts
    fprintf('\nNumber of clusters: \n Requested: %d; Returned: %d.\n', obj.num_clusts, this.k);
    obj.num_clusts = length(obj.clusters);
end


% Transpose the clusters
if obj.clustering_dim == 2 && obj.lpca_clustering_options == 0
    for i = 1 : obj.num_clusts
        obj.clusters{i} = obj.clusters{i}';
    end
end



%% Perform LPCA

% Warn that LPCA is being performed
fprintf('\nLPCA is being performed...\n\n');


% Create the options structure
lpca_options.cs = false;            
lpca_options.k = obj.k_init;   
lpca_options.con = obj.con;
lpca_options.islocal = true;
lpca_options.m = obj.m;         % Needed for CPCA
lpca_options.d = obj.d;         % Needed for CPCA


% Build the Local PCA objects
obj.lpca = {};
if ~isempty(obj.clusters)
    
    for i = 1 : obj.num_clusts
    
        % If you clustered on the first dimension, you need to pass the 
        % clustered mean and scaling factors, too
        if obj.clustering_dim == 1
            lpca_options.m = obj.m(obj.idx == i);
            lpca_options.d = obj.d(obj.idx == i);
        end
        
        % LPCA
        obj.lpca{i} = KPCA(obj.clusters{i}, lpca_options); % Create the KPCA objects
        obj.lpca{i}.k = obj.k; % You might use: obj.lpca{i}.getAppOrder(.9997);
        obj.lpca{i}.update(); % Trigger the updates
        obj.lpca{i}.trendFun = obj.trendFun; % Set the same trend function
        obj.lpca{i}.corrFun = obj.corrFun; % Set the same correlation function
        obj.lpca{i}.corrFun_optimization = obj.corrFun_optimization;
        
        % Useful for the Q-matrix
        if obj.clustering_dim == 2
            obj.lpca{i}.a_moresamples = obj.lpca{i}.modes' * obj.Ycs;
            obj.lpca{i}.moresamples = obj.xp;
        end
        
    end
    
end

fprintf('... LPCA was performed.');


end


function lpca = parallel_lpca(num_clusts, clustering_dim, m, d, ...
    clusters, k, trendFun, corrFun, corrFun_optimization, Ycs, xp)

lpca_options = cell(num_clusts, 1);
lpca = cell(num_clusts, 1);

parfor i = 1 : num_clusts
    % If you clustered on the first dimension, you need to pass the 
    % clustered mean and scaling factors, too
    if clustering_dim == 1
        lpca_options{i}.m = m(idx == i);
        lpca_options{i}.d = d(idx == i);
    end

    % LPCA
    lpca{i} = KPCA(clusters{i}, lpca_options{i}); % Create the KPCA objects
    lpca{i}.k = k; % You might use: obj.lpca{i}.getAppOrder(.9997);
    lpca{i}.update(); % Trigger the updates
    lpca{i}.trendFun = trendFun; % Set the same trend function
    lpca{i}.corrFun = corrFun; % Set the same correlation function
    lpca{i}.corrFun_optimization = corrFun_optimization;

    % Useful for the Q-matrix
    if clustering_dim == 2
        lpca{i}.a_moresamples = lpca{i}.modes' * Ycs;
        lpca{i}.moresamples = xp;
    end        
end

end





