function trainCpcaKriging(obj, varargin)

% Cannot run if no training points are present
if isempty(obj.training_points)
    warning('No training points found.');
    return
end

% Train the Kriging model for Global PCA
obj.cpcaKrigingModel = obj.train(obj.cpca_scores);

end
