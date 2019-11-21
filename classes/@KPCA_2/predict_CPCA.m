function [varargout] = predict_CPCA(obj, varargin)
% Prediction points
[n_obs, ~] = size(obj.prediction_points);
% Use Kriging predict()
n_scores = length(obj.cpcaKrigingModel); % Number of scores
obj.kriged_cpca_scores = zeros(n_scores, n_obs); % Initialize
for ii = 1 : n_scores 
    obj.kriged_cpca_scores(ii,:) = obj.cpcaKrigingModel{ii}.predict(obj.prediction_points);
end
% Get the predictions
obj.getKcpcaPredictions();
% Terminate here if this is a LPCA object
if obj.is_local
    return
end
% Estimate the errors
obj.getKcpcaPredictionsErrors();
end