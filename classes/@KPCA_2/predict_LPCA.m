function [varargout] = predict_LPCA(obj, varargin)
% Prediction points
[n_obs, ~] = size(obj.prediction_points);
% Use Kriging predict()
for ii = 1 : length(obj.local_pca)
    % All clusters will try to make the same predictions
    obj.local_pca{ii}.prediction_points = obj.prediction_points;
    % Predictions
    evalc('obj.local_pca{ii}.predict_PCA();');
end
% Get predictions of LPCA
obj.getKlpcaPredictions();
% Estimate the errors
obj.getKlpcaPredictionsErrors();
end



