function getKcpcaPredictions(obj, varargin)

% Get PCA approximation order
k = obj.pca_approximation_order;

% Check this is possible
if k > size(obj.pca_modes, 2)
    k = size(obj.pca_modes, 2);
end
if k > size(obj.kriged_cpca_scores, 1)
    k = size(obj.kriged_cpca_scores, 1);
end

% Decode scores and recover data(the pca_decode() method needs one
% additional input when decoding Local PCA scores)
if ~obj.is_local || (obj.is_local && obj.is_mesh_variable)
    obj.kcpca_predictions = obj.pca_decode(...
        obj.pca_modes(:,1:k), ...
        obj.kriged_cpca_scores(1:k,:), ...
        obj.scaling_factors, ...
        obj.mean_column);
else
    obj.kcpca_predictions = obj.pca_decode(...
        obj.pca_modes(:,1:k), ...
        obj.kriged_cpca_scores(1:k,:), ...
        obj.scaling_factors, ...
        obj.mean_column, ...
        obj.centroid);
end

end