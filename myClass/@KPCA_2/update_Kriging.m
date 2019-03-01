function [varargout] = update_Kriging(obj, first_input, second_input, varargin)
%% Pre-process
% Test data 
test_data = [];
if ~isempty(varargin)
    test_data = varargin{1};
end

% Default guess for the Kriging hyperparameters
guess_default = obj.hyp_guess;

% Samples and Values
samples = obj.training_points;              % Training points
values = first_input;                       % Training targets
prediction_points = obj.prediction_points;  % Prediction points

% If this Kriging is being applied on a LPCA (clustered observations)
if obj.is_local && obj.clustering_dimension && obj.includeMoreSamples
    samples = obj.moresamples;        % Training points
    values = obj.morevalues;          % Training targets
end

% Center-scale samples?
if obj.cs_samples
    [samples, mu, sf] = center_scale_data(samples', obj.scaling_criterion); % zscore(samples, 0, 1);
    samples = samples';
    if ~isempty(prediction_points)
        prediction_points = center_scale(prediction_points, mu, sf);
%         prediction_points = prediction_points';
    end
end

% Shorthen their names for convenience
trendFun = obj.trend_function;
corrFun = obj.correlation_function;

%% Kriging 
% Corr Fun Optimization
% if obj.correlation_function_optimization 
%     [I, guess] = corrFun_Optimize(samples, values, prediction_points, trendFun, obj.correlation_functions, test_data);
%     obj.correlation_function = I; % Set optimum correlation function
% % Then performKriging() will be re-run, with the optimum CORRFUN
% else
%     guess = guess_default;
% end

% In case of memory problems
[samples, values] = memory_test(samples, values, obj.is_mesh_variable, ...
    size(obj.mesh,1), size(obj.mesh,2), obj.selected_mesh_points);

% Priorities:
% [1] - Linear/MARS
% [2] - GPML
% [3] - Stacking
% [4] - Kriging

% Choose how to build the regression model
if obj.is_linear 
    % Polynomial interpolation (MARS)
    [regOutput, mse, k] = KPCA_2.linearRegression(samples, values,...
            prediction_points, obj.mars_maxOrd);
elseif obj.is_gpml
    % Use the GPML toolbox
    options.meanfunc = obj.meanfunc; % Leave empty for no mean function
    options.covfunc = obj.covfunc; % Covariance function
    options.likfunc = obj.likfunc; % Likelihood
    options.inffunc = obj.inffunc; % Inference method
    [regOutput, mse, k] = runGpml(samples, values', prediction_points, options, obj.hyp, obj.h_optimize);
    % The function runGpml() needs a transpose of VALUES.
else 
    if ~obj.stacking
        % Perform Kriging
        [regOutput, mse, k] = performKriging(samples, values,...
            prediction_points, trendFun, corrFun, obj.hyp_guess, obj.is_parallel);
    elseif iscell(obj.correlation_function)
        % Stacking on different Kriging models (user-supplied correlation
        % functions only)
        [regOutput, mse, k] = KPCA_2.ensembleKriging(samples, values,...
            prediction_points, obj.trend_functions, obj.correlation_function);
    else
        % Stacking on different Kriging models (all possible correlation
        % functions)
        [regOutput, mse, k] = KPCA_2.ensembleKriging(samples, values,...
            prediction_points, obj.trend_functions, obj.correlation_functions);
    end
end

obj.mse = mse;

%% Assign the output
% If there are no prediction points, the model will be trained but no
% prediction will be provided.
if nargout < 1
    if strcmp(second_input, 'pca') 
        obj.kriged_pca_scores = regOutput;
        obj.pcaKrigingModel = k;
        obj.pcaKrigingMSE = mse;
    elseif strcmp(second_input, 'direct')
        obj.kriged_direct_data = regOutput;
        obj.directKrigingModel = k;
        obj.directKrigingMSE = mse;
    elseif strcmp(second_input, 'cpca')
        obj.kriged_cpca_scores = regOutput;
        obj.cpcaKrigingModel = k;
        obj.cpcaKrigingMSE = mse;
    end
else
    varargout{1} = regOutput;
    if nargout > 1
        varargout{2} = k;
    end
    if nargout > 2
        varargout{3} = mse;
    end
end

end





