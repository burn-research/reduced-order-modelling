function runReconstruction(obj, varargin)
%% Notes:
% TO DO: include routines for the kriging models

%% Change app order
if ~isempty(varargin)
    obj.pca_approximation_order = varargin{1};
end

%% PCA
% Recover data and estimate errors
if ~isempty(obj.pca_scores)
    obj.recoverPca();
    obj.getPcaErrors();
end
if ~isempty(obj.cpca_scores)
    obj.recoverCpca();
    obj.getCpcaErrors();
end
if ~isempty(obj.local_pca) && ~isempty(obj.local_pca{1}.pca_scores)
    obj.recoverLpca();
    obj.getLpcaErrors();
end
if ~isempty(obj.local_pca) && ~isempty(obj.local_pca{1}.cpca_scores)
    obj.recoverLcpca();
    obj.getLcpcaErrors();
end

%% Kriging
% Recover predictions
if ~isempty(obj.kriged_pca_scores)
    obj.getKpcaPredictions();
end
if ~isempty(obj.kriged_cpca_scores)
    obj.getKcpcaPredictions();
end
if ~isempty(obj.local_pca) && ~isempty(obj.local_pca{1}.kriged_pca_scores)
    for ii = 1 : length(obj.local_pca)
        obj.local_pca{ii}.pca_approximation_order = obj.pca_approximation_order;
    end
    obj.getKlpcaPredictions();
end
if ~isempty(obj.local_pca) && ~isempty(obj.local_pca{1}.kriged_cpca_scores)
    for ii = 1 : length(obj.local_pca)
        obj.local_pca{ii}.pca_approximation_order = obj.pca_approximation_order;
    end
    obj.getKlcpcaPredictions();
end

% Estimate prediction errors
obj.getPredictionErrors();

end


