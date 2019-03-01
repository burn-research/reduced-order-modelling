function [mu, s2, hyp2] = Copy_of_runGpml(x, y, xs, varargin)
% NOTES:
% Trying to find the best way to control the guess h-parameters values. In
% this  function they are finally stored in the STRUCT 'hyp'.
% For now, if any difficulty arises, just do not provide the additional
% inputs to this function in the update_Kriging() routine, so that you can
% manually change the default values here.
%

n_args = length(varargin);

% Get useful sizes
y_dim = size(y,2);
n_points = size(xs,1);
x_dim = size(x,2);

% Setup important regression properties
options.meanfunc = {@meanSum, {@meanLinear, @meanConst}}; % Leave empty for no mean function
options.covfunc = @covSEard; % Covariance function
options.likfunc = @likGauss; % Likelihood
options.inffunc = @infGaussLik; % Inference method
if n_args > 0 && isstruct(varargin{1})
    options.meanfunc = varargin{1}.meanfunc; % Leave empty for no mean function
    options.covfunc = varargin{1}.covfunc; % Covariance function
    options.likfunc = varargin{1}.likfunc; % Likelihood
    options.inffunc = varargin{1}.inffunc; % Inference method
end

% Initialize the hyperparameters
hyp.mean = [ones(x_dim,1); 1];
hyp.cov = [0.1; ones(x_dim,1)];
hyp.lik = -1;
if n_args > 1
    hyp = varargin{2};
    if iscell(hyp) && length(hyp) ~= y_dim
        warning('[runGpml] You provided %d sets of h-parameters but there are %d scalars.', length(hyp), y_dim);
        warning('Trying to proceed using only hyp{1}.');
        hyp = hyp{1};
    end
end

% If more than one scalar output, use recursive function
if y_dim > 1
    % Initialize function outputs
    mu = zeros(n_points, y_dim);
    s2 = zeros(n_points, y_dim);
    hyp2 = cell(y_dim,1);
    % Re-run this same function, one scalar per run
    parfor ii = 1 : y_dim
        if iscell(hyp)
            [temp1, temp2, temp3] = Copy_of_runGpml(x, y(:,ii), xs, options, hyp{ii});
        else
            [temp1, temp2, temp3] = Copy_of_runGpml(x, y(:,ii), xs, options, hyp);
        end
        if n_points > 0
            mu(:,ii) = temp1;
            s2(:,ii) = temp2;
        end
        hyp2{ii} = temp3;
    end
    if n_points > 0
        mu = mu';
        s2 = s2';
    end
    % Do not forget to stop here
    return
end

% Get best set of hyperparameters
n_iter = 350;
evalc('hyp2 = minimize(hyp, @gp, -n_iter, options.inffunc, options.meanfunc, options.covfunc, options.likfunc, x, y);');

% Run the GP
if n_points > 0
    [mu, s2] = gp(hyp2, options.inffunc, options.meanfunc, options.covfunc, options.likfunc, x, y, xs);
else
    % If just for training, return empty prediction info. The hyp2
    % structure already has the optimized values of the hyper-paramters.
    mu = []; 
    s2 = [];
end

end

function [mu, s2, hyp2] = GPR(x, y, xs, options, hyp)
n_points = size(xs,1);
% Get best set of hyperparameters
n_iter = 350;
hyp2 = [];
evalc('hyp2 = minimize(hyp, @gp, -n_iter, options.inffunc, options.meanfunc, options.covfunc, options.likfunc, x, y);');

% Run the GP
if n_points > 0
    [mu, s2] = gp(hyp2, options.inffunc, options.meanfunc, options.covfunc, options.likfunc, x, y, xs);
else
    % If just for training, return empty prediction info. The hyp2
    % structure already has the optimized values of the hyper-paramters.
    mu = []; 
    s2 = [];
end
end



