function [krigingOutput, mse, k_models] = stackKriging(samples, values, p_points, t_funs, c_funs, guesses, varargin)
%% Description:
% More than one Kriging model will be trained, then they will be combined.
% Kriging models will differ not only in different trend and correlation
% functions (or set of hyperparameters) but may also differ in the training
% set.

%% Input processing
% Varargin{1}: Number of permutations or sets of samples
model_size = ceil(0.2 * size(samples,1));
n_models = nchoosek(size(samples,1), model_size); % Number of permutations
if ~isempty(varargin)
    if length(varargin{1}) == 1
        n_models = varargin{1};
    end
    if length(varargin) > 1 && length(varargin{2}) == 1
        model_size = varargin{2};
    end
end
% Check variable size and variable type: guesses
if size(guesses,1) ~= 1 && size(guesses,2) ~= 1
    error('The variable GUESSES must be a vector or a scalar.');
else
    n_guess = length(guesses);
end
% Check variable size and variable type: c_funs
if isa(c_funs, 'function_handle')
    n_corr = 1;
    temp = {c_funs}; c_funs = temp; % Convert it to cell-array for convenience
elseif isa(c_funs, 'cell')
    n_corr = length(c_funs);
else
    error('The variable C_FUNS must be a function handle or a cell array.');
end
% Check variable size and variable type: t_funs
if isa(t_funs, 'char')
    n_trend = 1;
    temp = {t_funs}; t_funs = temp; % Convert it to cell-array for convenience
elseif isa(t_funs, 'cell')
    n_trend = length(t_funs);
else
    error('The variable T_FUNS must be a function handle or a cell array.');
end
% Predictions size
[n_var, n_samples] = size(values); % values: (variables x samples)

%% Get sets for the Kriging models
[setOfSamples, setOfValues, setOfPoints, setOfTargets] = ...
    get_permSets(n_models, model_size, samples, values);

%% Train a superior Kriging model for each variable
n_points = size(p_points, 1);
krigingOutput = zeros(n_var, n_points);
for j = 1 : n_var
    % Create n_models inferior Kriging models for variable j
    k_models = cell(n_models, 1);
    for ii = 1 : n_models
        [~, ~, k_models{ii}] = KPCA_2.performKriging(setOfSamples{ii}, setOfValues{ii}(j,:), [], ...
            t_funs{1}, c_funs{1}, guesses(1));
    end

    % For each setOfPoints{ii}, get the models' predictions'
    models_predictions = zeros(n_models, n_samples);
    for ii = 1 : n_models
        models_predictions(ii,:) = k_models{ii}.predict(samples);
    end
    targets = values(j,:);
    
    % Train one superior Kriging model for variable j
    [~, ~, v_model] = KPCA_2.performKriging(models_predictions', targets, [], ...
        t_funs{1}, c_funs{1}, guesses(1));
    
    % Get predictions of each small model for variable j on p_points
    y_out = zeros(n_points, n_models);
    for ii = 1 : n_models
        y_out(:,ii) = k_models{ii}.predict(p_points);
    end
    
    krigingOutput(j,:) = v_model.predict(y_out);
end

% Output
mse = [];
k_models = [];

end


function [setOfSamples, setOfValues, setOfPoints, setOfTargets] = get_permSets(n_models, model_size, samples, values)
% NOTES:
% [~, idx] = datasample(1:100, 4, 'Replace', false);
%

% Sizes
n_samples = size(samples, 1);

% 
idx = cell(n_models,1);
setOfSamples = cell(n_models,1);
setOfValues = cell(n_models,1);
setOfPoints = cell(n_models,1);
setOfTargets = cell(n_models,1);
if n_models == 1
    setOfSamples{1} = samples;
    setOfValues{1} = values;
    return
end

% Check number of sought permutations
temp_vec = round(1:n_samples);
for ii = 1 : n_models
    [setOfPoints{ii}, idx{ii}] = datasample(samples, model_size, 1, 'Replace', false);
    setOfTargets{ii} = values(:,idx{ii});
    I = ~ismember(temp_vec, idx{ii});
    setOfSamples{ii} = samples(I,:);
    setOfValues{ii} = values(:,I);
end

end



