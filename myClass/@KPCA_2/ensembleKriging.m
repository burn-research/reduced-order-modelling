function [krigingOutput, mse_out, k_models] = ensembleKriging(samples, values, prediction_points, trendFuns, corrFuns, varargin)
%% Description:
% More than one Kriging model will be trained, then they will be combined.
% Kriging models will differ only in different trend and correlation
% functions. They will not differ in the training set nor in the
% hyperparameters.

%% Input processing
[nt, nc, dim_in, dim_out, trendFuns, corrFuns] = ...
    process_input(samples, values, prediction_points, trendFuns, corrFuns);

%% Train models
n_models = nt * nc; % Number of models to train
kk = 0; % Iterator
predictions = cell(n_models,1); % Predicitons of each model
mse = cell(n_models,1); % MSE of each model
k_models = cell(n_models,1); % Cell-array to store every trained model
mse_model = zeros(n_models,1); % Measures of each model performance
% Train a model for each correlation fun
for ii = 1 : nc
    % Train a model for each trend fun
    for jj = 1 : nt
        [~, guess] = corrFun_Optimize(samples, values, prediction_points, ...
            trendFuns{jj}, corrFuns{ii});
        kk = kk + 1; % Increment the iterator
        try
            % Perform Kriging
            [predictions{kk}, mse{kk}, k_models{kk}] = performKriging(samples, values, ...
                prediction_points, trendFuns{jj}, corrFuns{ii}, guess);
            % For each model, get a general MSE
            temp = abs(mse{kk}) ./ (abs(predictions{kk}) + eps);
            mse_model(kk) = max( mean( temp, 1) );
        catch ME 
            % In case the algorithm fails, return bad outputs
            n = size(values,1); m = size(prediction_points,1);
            predictions{kk} = zeros(n,m); 
            mse{kk} = zeros(n,m) + Inf;
            % For each model, get a general MSE
            mse_model(kk) = Inf;
        end
    end
end
% Weights inversely proportional to MSE
w = 1 ./ mse_model; 
w = w ./ sum(w);

%% Output
% If no prediction points are present, just return the kriging models
if isempty(prediction_points)
    return
end
% Now we need to combine the different models' predictions
n_points = size(prediction_points, 1); % Number of prediction points
krigingOutput = zeros(dim_out, n_points);
% Return a weighted MSE using the weights 
mse_out = zeros(dim_out, n_points);
for ii = 1 : n_models
    if mse{ii} ~= Inf
        mse_out = mse_out + mse{ii} * w(ii);
    else
        mse_out = mse_out + 0;
    end
end
% For each prediction point, combine each model's prediction for that point
for ii = 1 : n_points
    % Get the weights only for this prediction point
    mse_model = zeros(n_models, 1); % Measures of each model performance
    for kk = 1 : n_models
        % MSE check
        temp = abs(mse{kk}(:,ii)) ./ (abs(predictions{kk}(:,ii)) + eps);
        mse_model(kk) = max(temp);
    end
    % Weights inversely proportional to MSE
    w = 1 ./ mse_model; 
    w = w ./ sum(w);
    % Use the weights to combine predictions
    sol = zeros(dim_out, n_models); % Columns are each model's prediction for one prediction point 
    for kk = 1 : n_models
        sol(:,kk) = predictions{kk}(:,ii); % Load prediction for model kk
    end
    krigingOutput(:,ii) = sol * w; % Weighed average of the predictions
end

end


function [nt, nc, dim_in, dim_out, trendFuns, corrFuns] = process_input(samples, values, prediction_points, trendFuns, corrFuns, varargin)
% Hiding the input processing here for readability

% Get input and output space dimensions and check the number of samples is
% equal to the number of observations
[n_samples, dim_in] = size(samples);
[dim_out, temp] = size(values);
if temp ~= n_samples
    error('First input`s number of rows must be equal to second input`s number of columns.');
end

% If present, do the check for prediction_points too
if ~isempty(prediction_points) && size(prediction_points,2) ~= size(samples,2)
    error('Training and prediction points must have the same dimension.');
end

% Let's check how many models we need to train
% Number of trend functions
if iscell(trendFuns)
    nt = length(trendFuns);
else
    trendFuns = {trendFuns}; % Make it a cell for uniformity
    nt = 1; 
end
% Number of correlation functions
if iscell(corrFuns)
    nc = length(corrFuns);
else
    corrFuns = {corrFuns}; % Make it a cell for uniformity
    nc = 1; 
end

end





