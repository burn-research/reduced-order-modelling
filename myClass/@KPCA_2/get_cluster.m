function [cluster, varargout] = get_cluster(obj, cluster_index)
%% Description
%{
    Retrieve a cluster by inputting its index. Output is a (R X M) matrix,
    with R being the number of rows (elements) belonging to the cluster, M
    being the number of columns of the original data matrix.
%}

%% Main
I = obj.local_nz_idx_clust{cluster_index}; % Get the needed rows
cluster = obj.training_data(I,:); % Extract those rows

end


