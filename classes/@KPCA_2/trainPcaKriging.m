function trainPcaKriging(obj, varargin)

% Cannot run if no training points are present
if isempty(obj.training_points)
    warning('No training points found.');
    return
end

% Train the Kriging model for Global PCA
obj.pcaKrigingModel = obj.train(obj.pca_scores);

end
