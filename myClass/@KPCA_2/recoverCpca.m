function recoverCpca(obj, varargin)

% Get PCA approximation order
k = obj.pca_approximation_order;

% Check this is possible
if k > size(obj.pca_modes, 2)
    k = size(obj.pca_modes, 2);
end
if k > size(obj.cpca_scores, 1)
    k = size(obj.cpca_scores, 1);
end

% Recover data
obj.cpca_recovered_data = obj.pca_decode(...
    obj.pca_modes(:,1:k), ...
    obj.cpca_scores(1:k,:), ...
    obj.scaling_factors, ...
    obj.mean_column...
    );

end

