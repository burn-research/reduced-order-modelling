function trainLcpcaKriging(obj, varargin)

% Cannot run if no training points are present
if isempty(obj.training_points)
    warning('No training points found.');
    return
end

% Cannot run if LPCA has not been performed yet
if isempty(obj.local_pca)
    warning('No LPCA found.');
    return
end

% Train the Kriging model for Local PCA
for ii = 1 : length(obj.local_pca)
    obj.local_pca{ii}.trainCpcaKriging();
end

end

