function [varargout] = predict_PCA(obj, varargin)
% Prediction points
[n_obs, ~] = size(obj.prediction_points);
% Use Kriging predict()
n_scores = length(obj.pcaKrigingModel); % Number of scores
obj.kriged_pca_scores = zeros(n_scores, n_obs); % Initialize
for ii = 1 : n_scores 
    obj.kriged_pca_scores(ii,:) = obj.pcaKrigingModel{ii}.predict(obj.prediction_points);
end
% Get the predictions
obj.getKpcaPredictions();
% Terminate here if this is a LPCA object
if obj.is_local
    return
end
% Estimate the errors
obj.getKpcaPredictionsErrors();
end



