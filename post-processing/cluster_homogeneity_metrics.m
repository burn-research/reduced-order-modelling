function [mean_nrmse, mean_silhouette, mean_delta, var_delta] = cluster_homogeneity_metrics(X, idx, q, cent_crit, scal_crit)
% This function computes four metrics of cluster homogeneity.
%
% Input:
% ------------
% - X
%       raw data matrix, uncentered, unscaled.
%
% - idx
%       cluters indices.
%
% - q
%       number of PCs to retain in Local PCA.
%
% - cent_crit
%       centering criteria as per `center` function.
%
% - scal_crit
%       scaling criteria as per `scale` function.
%
% Output:
% ------------
% - mean_nrmse
%       mean normalized Root Mean Squared Error (RMSE).
%
% - mean_silhouette
%       mean silhoutte.
%
% - mean_delta
%       mean of dy/mean(y).
%
% - var_delta
%       variance of dy/mean(y), where y are the observations belonging to a particular cluster.

%% cluster_homogeneity_metrics()
% Center and scale data:
[X_cent_scal, centerings] = center(X, cent_crit);
[X_cent_scal, scalings] = scale(X_cent_scal, X, scal_crit);

% Perform PCA in local cluster and recover the data set:
[clusters] = get_clusters(X_cent_scal, idx);
[eigvec, scores, eigenvalues, centroids, local_scalings] = lpca(clusters, q, cent_crit);
[X_rec] = recover_from_lpca(idx, eigvec, scores, q, centroids, local_scalings, centerings, scalings);

% Mean normalized Root Mean Squared error:
[~, nrmse] = quality_of_reconstruction_measures(X, X_rec);
mean_nrmse = mean(nrmse);

% Mean silhoutte:
mean_silhouette = mean(silhouette(X, idx));

% dy/mean(y) factor, mean and variance:
n_var = size(X, 2);
n_clust = size(clusters, 1);

delta = zeros(n_clust, n_var);

groups = get_clusters(X, idx);

for i = 1:1:n_clust
    for j = 1:1:n_var
        delta(i,j) = (max(groups{i}(:,j)) - min(groups{i}(:,j)))/(mean(groups{i}(:,j) + 1e-16));
    end
end

tmp = zeros(1,n_clust);

for i = 1:1:n_clust
    tmp(i) = mean(delta(i,:));
end

mean_delta = mean(tmp);
var_delta = var(tmp);

end
