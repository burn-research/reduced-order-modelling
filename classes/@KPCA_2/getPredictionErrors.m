function getPredictionErrors(obj, varargin)

% PCA
if ~isempty(obj.kpca_predictions)
    obj.getKpcaPredictionsErrors();
end
% CPCA
if ~isempty(obj.kriged_cpca_scores)
    obj.getKcpcaPredictionsErrors();
end
% LPCA
if ~isempty(obj.local_pca) && ~isempty(obj.local_pca{1}.kriged_pca_scores)
    obj.getKlpcaPredictionsErrors();
end
% LCPCA
if ~isempty(obj.local_pca) && ~isempty(obj.local_pca{1}.kriged_cpca_scores)
    obj.getKlcpcaPredictionsErrors();
end
% Direct
if ~isempty(obj.kriged_direct_data)
    obj.getDirectKrigingErrors();
end

end