function [combined_predictions, varargout] = combinePredictions(predictions, mse, kriged, varargin)

if length(predictions) ~= length(mse)
    error('Number of models should be the same.');
elseif length(mse) ~= length(kriged)
    error('Number of models should be the same.');
end

% User can provide distances
dist = [];
if ~isempty(varargin)
    dist = varargin{1};
end

% Get useful info
[dim_out, n_predictions] = size(predictions{1});
n_models = length(predictions);

% Rebuild for each prediction point
combined_predictions = zeros(dim_out, n_predictions);
for jj = 1 : n_predictions
    % For each cluster, get the j-th prediction and weigh it
    % accordingly. The prediction for the j-th prediction point will be
    % a weighed average of the clusters predictions, with the distance
    % evaluated before as weights
    % For prediction jj, get one MSE per model/cluster
    mse_model = zeros(n_models, 1);
    for ii = 1 : n_models
        % MSE check
        temp = abs(mse{ii}(:,jj)) ./ (kriged{ii}(:,jj) + eps);
        mse_model(ii) = max( mean( temp, 1) ); % MSE for the ii-th model
    end
    % Weights inversely proportional to MSE
    weights = 1 ./ mse_model; 
    weights = weights ./ sum(weights); % Weights for the jj-th prediction
    sol = zeros(dim_out, n_models); % Columns are each model's prediction for jj-th prediction point 
    for ii = 1 : n_models
        sol(:,ii) = predictions{ii}(:,jj); % Load prediction for model kk
    end
    combined_predictions(:,jj) = sol * (weights); % Weighed average of the predictions
    % Use also the distances as meta-features
    if ~isempty(dist)
        weights_2 = 1 ./ dist(:,jj);
        weights_2 = weights_2 ./ sum(weights_2);
        combined_predictions(:,jj) = sol * (weights + weights_2) / 2;
    end
    % We can also get the best one, instead of combining
%         [~, I] = max(weights);
%         obj.klpca_predictions(:,jj) = sol(:,I);
end

end


