function localPCA_accelerated(obj, varargin)
%% Description
% Accelerated clustering: maximun number of levels is 2. The user asked for
% a number of clusters equal to (obj.number_of_clusters). Data is divided
% into sqrt(obj.number_of_clusters) clusters, and then each cluster is again
% partinioned into sqrt(obj.number_of_clusters) clusters. 
% The number of levels is kept 2 for convenience, as it might already
% provide enough acceleration.
%

%% Accelerated Clustering
% According to how this routine is called, this routine has to return the
% definitive obj.local_idx that will be then used when calling the routine
% obj.runLpca().

% Requested number of clusters
requested_n_clusters = obj.number_of_clusters;

% The number of levels is 2
n_levels = 2;

% Cell-array where to store the number of clusters for each level
b = cell(n_levels, 1);

% Number of clusters at the first level
b{1} = ceil(sqrt(obj.number_of_clusters));

% Run clustering on the first level
if ~obj.clustering_dimension
    [obj.local_idx, ~, this] = localPCA(obj.centered_scaled_data, obj.pca_approximation_order, b{1}, obj.idx_0);
else
    [obj.local_idx, ~, this] = localPCA(obj.centered_scaled_data', obj.pca_approximation_order, b{1}, obj.idx_0);
end
obj.local_nz_idx_clust = this.nz_idx_clust;

% Check the number of clusters returned is the same
if length(obj.local_nz_idx_clust) ~= b{1}
    % Update the number of clusters
    obj.number_of_clusters = length(obj.local_nz_idx_clust);
    b{1} = length(obj.local_nz_idx_clust);
end

% Initialize the LPCA objects for the first level
obj.local_pca = cell(b{1}, 1);

% Number of clusters at the second level
b{2} = ceil(obj.number_of_clusters / b{1});

% Run clustering on the second level
clusters_count = 0;
for ii = 1 : b{1}
    % Create options
    options.mean_column = obj.mean_column(obj.local_nz_idx_clust{ii});
    options.scaling_factors = obj.scaling_factors(obj.local_nz_idx_clust{ii});
    options.is_local = true;
    
    % Create Local PCA objects
    obj.local_pca{ii} = KPCA_2( obj.centered_scaled_data(obj.local_nz_idx_clust{ii}, :), options );
    
    % Run the clustering inside the LPCA objects
    obj.local_pca{ii}.number_of_clusters = b{2};
    evalc('obj.local_pca{ii}.runClustering();');
    
    % Update the final local_idx: offset for the cluster index. 'local_idx'
    % is updated by 'local_pca{ii}.local_idx', every Local PCA object looks
    % for the same number of clusters: their indeces' range is the same. To
    % differentiate them, the first set of indices does not have an offset
    % (clusters_count = 0 before the for loop), then it increases so that,
    % for example, for the second for-loop iteration, the lowest index
    % value for 'local_pca{2}.local_idx' is 1 +
    % local_pca{1}.number_of_clusters.
    obj.local_idx(obj.local_nz_idx_clust{ii}) = obj.local_pca{ii}.local_idx + clusters_count; 
    
    % Update the clusters count
    clusters_count = clusters_count + obj.local_pca{ii}.number_of_clusters;
end

% Update the final local_nz_idx_clust
obj.local_nz_idx_clust = here_partitionVQ(obj.local_idx);

% Check the number of clusters returned is the same
if length(obj.local_nz_idx_clust) ~= requested_n_clusters
    fprintf('\nNumber of clusters: \n Requested: %d; Returned: %d.\n', ...
        requested_n_clusters, length(obj.local_nz_idx_clust));
end

% Update the number of clusters
obj.number_of_clusters = length(obj.local_nz_idx_clust);

% Clear the temporary local PCA objects
obj.local_pca = {};

if length(obj.local_nz_idx_clust) ~= length(unique(obj.local_idx))
    error('Something went wrong!');
end


end


function [nz_idx_clust, varargout] = here_partitionVQ(idx, varargin)

    k = max(idx);

    idx_clust = cell(k, 1);
    for j = 1 : k
        idx_clust{j} = find(idx == j);
    end
    
    n_points = zeros(k, 1); 
    for j = 1 : k
        idx_clust{j} = find(idx == j);
        n_points(j) = size(idx_clust{j}, 1);
    end
    nz_idx = find(n_points > -1);
    
    nz_idx_clust = cell(k, 1);
    for j = 1 : k
        nz_idx_clust{j} = idx_clust{nz_idx(j)};
    end
    
end



