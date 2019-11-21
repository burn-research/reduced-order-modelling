function runClustering(obj, varargin)
%% Description
% Routine that partitions data into clusters.
%

%% Input
n_args = length(varargin);
if n_args > 0
    obj.number_of_clusters = varargin{1};
end

%% Check the number of clusters
if obj.number_of_clusters == 1
    obj.local_idx = ones(size(obj.training_data, obj.clustering_dimension + 1), 1);
    obj.local_pca = {};
    obj.local_nz_idx_clust = {};
    fprintf('\nLocal PCA was run, with only one cluster.\n');
    return; % Do not run the clustering optimization procedure
end

%% Run the appropriate clustering routine
fprintf('\nClustering data... \n');
if ~obj.accelerated_clustering
    % Choose if to cluster variables or samples
    if ~obj.clustering_dimension 
        % Cluster variables
        obj.local_idx = localPCA(...
            obj.centered_scaled_data, ...
            obj.pca_approximation_order, ...
            obj.number_of_clusters, ...
            obj.idx_0, ... % Centroids initialization
            obj.lpca_clustering_options); 
        % Check the number of clusters returned is the same
        returned_k = length(unique(obj.local_idx));
        if returned_k ~= obj.number_of_clusters
            fprintf('\nNumber of clusters: \n Requested: %d; Returned: %d.\n', ...
                obj.number_of_clusters, returned_k);
        end
        % Update the number of clusters
        obj.number_of_clusters = returned_k; 
    else
        % Cluster observations
        obj.local_idx = localPCA(...
            obj.centered_scaled_data', ...
            obj.pca_approximation_order, ...
            obj.number_of_clusters, ...
            obj.idx_0, ... % Centroids initialization
            obj.lpca_clustering_options);
        % Check the number of clusters returned is the same
        returned_k = length(unique(obj.local_idx));
        if returned_k ~= obj.number_of_clusters
            fprintf('\nNumber of clusters: \n Requested: %d; Returned: %d.\n', ...
                obj.number_of_clusters, returned_k);
        end
        % Update the number of clusters
        obj.number_of_clusters = returned_k;
    end   
else
    obj.localPCA_accelerated();
end
fprintf('\n...data clustered. \n');
% Remove holes in IDX
obj.local_idx = removeHoles(obj.local_idx);
% Save
obj.local_idx_stored{end+1} = obj.local_idx;
obj.number_of_clusters_stored(end+1) = numel(unique(obj.local_idx));

end




