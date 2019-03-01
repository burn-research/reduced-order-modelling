function choosePartitioningCriteria(obj, varargin)
%% Run the appropriate clustering routine
if ~obj.accelerated_clustering
    % Choose if to cluster variables or samples
    if ~obj.clustering_dimension 
        % Cluster variables
        [obj.local_idx, ~, this] = obj.localPCA(...
            obj.centered_scaled_data, ...
            obj.pca_approximation_order, ...
            obj.number_of_clusters, ...
            obj.idx_0); % Centroids initialization
        % Local clusters
        obj.local_nz_idx_clust = this.nz_idx_clust;
        % Check the number of clusters returned is the same
        if length(obj.local_nz_idx_clust) ~= obj.number_of_clusters
            fprintf('\nNumber of clusters: \n Requested: %d; Returned: %d.\n', ...
                obj.number_of_clusters, this.k);
        end
        % Update the number of clusters
        obj.number_of_clusters = length(obj.local_nz_idx_clust); 
    else
        % Cluster observations
        [obj.local_idx, ~, this] = obj.localPCA(...
            obj.centered_scaled_data', ...
            obj.pca_approximation_order, ...
            obj.number_of_clusters, ...
            obj.idx_0);
        % Local clusters
        obj.local_nz_idx_clust = this.nz_idx_clust;
        % Check the number of clusters returned is the same
        if length(obj.local_nz_idx_clust) ~= obj.number_of_clusters
            fprintf('\nNumber of clusters: \n Requested: %d; Returned: %d.\n', ...
                obj.number_of_clusters, this.k);
        end
        % Update the number of clusters
        obj.number_of_clusters = length(obj.local_nz_idx_clust);
    end   
else
    obj.localPCA_accelerated();
end

end


