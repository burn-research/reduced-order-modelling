function [predictions, mse, k] = performKriging(X, Y, Xs, trendFun, corrFun, guess, is_parallel, varargin)
%% Description
% X  (obs x var): training points
% Y  (var x obs): training values
% Xs (obs x var): [optional] prediction points
%

%% Input
% Trend function
if ~exist('trendFun', 'var') || isempty(trendFun)
    trendFun = 'regpoly0';
end
% Correlation function
if ~exist('corrFun', 'var') || isempty(corrFun)
    corrFun = @corrgauss;
end
% Starting guess for the evaluation of the hyperparameters for Kriging
if ~exist('guess', 'var') || isempty(guess)
    guess = [];
end
% Parallel
if ~exist('is_parallel', 'var') || isempty(is_parallel)
    is_parallel = true;
elseif ~islogical(is_parallel)
    is_parallel = (is_parallel > 0);
end

%% Main
% Generate starting guess for the hyper-parameters
if isempty(guess)
    fprintf('Generating hyper-parameters...\n');
    guess = generateHyperparameters0(X, Y, trendFun, corrFun, false);
    fprintf('Hyper-parameters generated.\n');
end
% We need to perform Kriging one scalar at time.
y_dim = size(Y, 1);
n_points = size(Xs, 1);
if n_points > 0
    predictions = zeros(y_dim, n_points);
    mse = zeros(y_dim, n_points);
    k = cell(y_dim, 1);
    if is_parallel
        parfor ii = 1 : y_dim
            if ~iscell(guess)
                [m, s, kr] = ood_level1(X, Y(ii,:), ...
                    Xs, trendFun, corrFun, guess);
            else
                [m, s, kr] = ood_level1(X, Y(ii,:), ...
                    Xs, trendFun, corrFun, guess{ii});
            end
            predictions(ii,:) = m;
            mse(ii,:) = s;
            k{ii} = kr;
        end
    else
        for ii = 1 : y_dim
            if ~iscell(guess)
                [m, s, kr] = ood_level1(X, Y(ii,:), ...
                    Xs, trendFun, corrFun, guess);
            else
                [m, s, kr] = ood_level1(X, Y(ii,:), ...
                    Xs, trendFun, corrFun, guess{ii});
            end
            predictions(ii,:) = m;
            mse(ii,:) = s;
            k{ii} = kr;
        end
    end
else
    predictions = [];
    mse = [];
    k = cell(y_dim, 1);
    if is_parallel
        parfor ii = 1 : y_dim
            if ~iscell(guess)
                [~, ~, kr] = ood_level1(X, Y(ii,:), ...
                    Xs, trendFun, corrFun, guess);
            else
                [~, ~, kr] = ood_level1(X, Y(ii,:), ...
                    Xs, trendFun, corrFun, guess{ii});
            end
            k{ii} = kr;
        end
    else
        for ii = 1 : y_dim
            if ~iscell(guess)
                [~, ~, kr] = ood_level1(X, Y(ii,:), ...
                    Xs, trendFun, corrFun, guess);
            else
                [~, ~, kr] = ood_level1(X, Y(ii,:), ...
                    Xs, trendFun, corrFun, guess{ii});
            end
            k{ii} = kr;
        end
    end
end

end


function [predictions, mse, k] = ood_level1(samples, values, prediction_points, trendFun, corrFun, guess)
% Kriging, training the model
k = ood(samples, values, guess, trendFun, corrFun);
k = k.cleanup();
if ~isempty(prediction_points)
    % Predictions
    [predictions, mse] = k.predict(prediction_points);
    % NaN check
    predictions = predictions';     mse = mse';
    checkNaN = isnan(predictions);  predictions(checkNaN) = 0;
    checkNaN = isnan(mse);          mse(checkNaN) = 0;
else
    predictions = [];
    mse = [];
end
end


function k = ood(samples, values, guess, varargin)
%% Description:
% Remember that this function is now onl accessed by scalar targets: the
% variables VALUES will contain observations of one scalar target.

%% Inputs
% Trend and correlation functions start out empty: then the user can
% provide an istance of them. 
trendFun = [];
corrFun = [];
n_args = length(varargin);
% User-provided trend function
if n_args > 0 
    trendFun = varargin{1}; 
end
% User-provided correlation function
if n_args > 1
    corrFun = varargin{2};  
end
% Default functions: if they are still empty (no user-provided value or
% empty value), they get assigned some default values
if isempty(trendFun)
    trendFun = 'regpoly2'; % Trend function: can be downgraded later on automatically
end
if isempty(corrFun)
    corrFun = @corrmatern32; % Correlation function
end

%% Setting important properties
% Check validity of the inputs
if ~strcmp(trendFun(1:end-1),'regpoly') && ~strcmp(trendFun,'')
    fprintf('\tRegression function was mispecified, the default one was set.\t');
    trendFun = 'regpoly2'; % If the trend fun was not properply specified, choose the default one
end
% Samples
values = values'; % Transpose VALUES (to be interfaced with the Kriging toolbox)
[nSamples, dim] = size(samples);
% Check if the chosen trend function is ok (should do the check for
% regpoly3, regpoly4 ecc)
if strcmp(trendFun, 'regpoly2')
    p = .5 * (dim + 1) * (dim + 2);
    if nSamples < p
        trendFun = 'regpoly1';
    end
end
if strcmp(trendFun, 'regpoly1')
    p = dim + 1;
    if nSamples < p
        trendFun = 'regpoly0';
    end
end
% Corrgaussp needs one more hyperparameter
if strcmp(corrFun, 'corrgaussp')
   dim = dim + 1; 
end

%% Kriging
% Kriging options
opts = Kriging.getDefaultOptions();
opts.generateHyperparameters0 = false;
opts.hpOptimizer = SQPLabOptimizer(dim, 1);
opts.regressionMaxOrder = dim; 
% Hyperparameters 
Val = 1e10;                            
lb = -Val * ones(1, dim);  ub = Val * ones(1, dim); % (1 x d) Lower and upper bounds
l_guess = length(guess);
if l_guess == dim
    theta0 = guess;
elseif l_guess == 1
    theta0 = guess * ones(1, dim);  % (1 x d) Guess
else
    % The user gave a wrong dimension to GUESS, we might then either
    % take the first (guess(1)) or mean(guess)
    theta0 = guess(1) * ones(1, dim);  % (1 x d) Guess
end
opts.hpBounds = [lb; ub]; % (2 x d) Hyperparameter optimization bounds
% Kriging
try
    evalc('k = Kriging(opts, theta0, trendFun, corrFun);');
    evalc('k = k.fit(samples, values); k = k.cleanup();');
catch ME 
    disp(ME);
    try 
        clear k; evalc('k = Kriging(opts, theta0, trendFun, corrFun);');
        evalc('k = k.fit(samples, values); k = k.cleanup();');
    catch ME 
        disp(ME);
        % If still an error occurs, just find a way to return a Kriging object k, even if meaningless
        clear k; evalc('k = Kriging(opts, theta0, trendFun, corrFun);');
        evalc('k = k.fit(samples, values * 0); k = k.cleanup();');
    end 
end

end

%     % Stochastic Kriging options
%     opts = Stochastic.getDefaultOptions();
%     opts.hpOptimizer = SQPLabOptimizer( dim, 1 );
%     opts.regressionMaxOrder = dim; 
%     % Hyperparameters 
%     Val = 1e10;                            
%     lb = zeros(1, dim);  ub = Val * ones(1, dim); % (1 x d) Lower and upper bounds
%     l_guess = length(guess);
%     if l_guess == dim
%         theta0 = guess;
%     elseif l_guess == 1
%         theta0 = guess * ones(1, dim);  % (1 x d) Guess
%     else
%         % The user gave a wrong dimension to GUESS, we might then either
%         % take the first (guess(1)) or mean(guess)
%         theta0 = guess(1) * ones(1, dim);  % (1 x d) Guess
%     end
%     opts.hpBounds = [lb; ub]; % (2 x d) Hyperparameter optimization bounds
%     % Blind Kriging
%     evalc('k = BlindKriging( opts, theta0, trendFun, corrFun );');
%     try
%         evalc('k = k.fit(samples, values);');
%     catch ME 
%         evalc('k = k.fit(samples, values);');
%     end
