function [mu, s2, hyp2, fmu, fs2, lp, post] = runGpml(x, y, xs, varargin)
%% Description
% Trying to find the best way to control the guess h-parameters values. In
% this function they are finally stored in the STRUCT 'hyp'.
% For now, if any difficulty arises, just do not provide the additional
% inputs to this function in the update_Kriging() routine, so that you can
% manually change the default values here.
% Inputs:
% - x, n_obs x n_var
% - y, n_obs x n_var
%

%% Input
n_args = length(varargin);
in_count = 0;
% Get useful sizes
y_dim = size(y,2);
n_points = size(xs,1);
x_dim = size(x,2);
% Setup important regression properties
options.meanfunc = {@meanSum, {@meanLinear, @meanConst}}; % Leave empty for no mean function
options.covfunc = @covSEard; % Covariance function
options.likfunc = @likGauss; % Likelihood
options.inffunc = @infGaussLik; % Inference method
if n_args > in_count && isstruct(varargin{in_count+1})
    in_count = in_count + 1;
    options.meanfunc = varargin{in_count}.meanfunc; % Leave empty for no mean function
    options.covfunc = varargin{in_count}.covfunc; % Covariance function
    options.likfunc = varargin{in_count}.likfunc; % Likelihood
    options.inffunc = varargin{in_count}.inffunc; % Inference method
end
% Initialize the hyperparameters
hyp.mean = ones(x_dim+1, 1);
hyp.cov = [0.1; ones(x_dim,1)];
hyp.lik = -1;
if n_args > in_count && ~isempty(varargin{in_count+1})
    in_count = in_count + 1;
    hyp = varargin{in_count};
    if iscell(hyp) && length(hyp) ~= y_dim
        warning('[runGpml] You provided %i sets of h-parameters but there are %i scalars.', length(hyp), y_dim);
        warning('Trying to proceed using only hyp{1}.');
        hyp = hyp{1};
    end
end
% Optimize hyperparamenters
if n_args > in_count
    in_count = in_count + 1;
    h_optimize = varargin{in_count};
    if h_optimize
        fprintf('Generating hyper-parameters...\n');
        [ls, beta] = generateHyperparameters0(x, y', options, hyp, false, 'gpr');
        fprintf('Hyper-parameters generated.\n');
        hyp_out = cell(length(ls),1);
        for ii = 1 : length(hyp_out)
            hyp_out{ii}.mean = beta{ii};
            hyp_out{ii}.cov = ls{ii};
            if iscell(hyp)
                hyp_out{ii}.lik = hyp{ii}.lik;
            else
                hyp_out{ii}.lik = hyp.lik;
            end
        end
        hyp = hyp_out; clear hyp_out;
    end
end
% Prediction only?
pred_only = false;
if n_args > in_count
    in_count = in_count + 1;
    pred_only = varargin{in_count};
end
%% Main
% Initialize function outputs
mu = zeros(n_points, y_dim);
s2 = zeros(n_points, y_dim);
hyp2 = cell(y_dim,1);
fmu = cell(y_dim,1);
fs2 = cell(y_dim,1);
lp = cell(y_dim,1);
post = cell(y_dim,1);
% One scalar per run
for ii = 1 : y_dim
    if iscell(hyp)
        [mu_ii, s2_ii, hyp2_ii, fmu_ii, fs2_ii, lp_ii, post_ii] = GPR(x, y(:,ii), xs, options, hyp{ii}, pred_only);
    else
        [mu_ii, s2_ii, hyp2_ii, fmu_ii, fs2_ii, lp_ii, post_ii] = GPR(x, y(:,ii), xs, options, hyp, pred_only);
    end
    % If n_points == 0, mu and s2 are empty
    if n_points > 0
        mu(:,ii) = mu_ii;
        s2(:,ii) = s2_ii;
    end
    hyp2{ii} = hyp2_ii;
    fmu{ii} = fmu_ii;
    fs2{ii} = fs2_ii;
    lp{ii} = lp_ii;
    post{ii} = post_ii;
end
% If n_points == 0, mu and s2 are empty
% if n_points > 0
%     mu = mu';
%     s2 = s2';
% end


end

function [mu, s2, hyp2, fmu, fs2, lp, post] = GPR(x, y, xs, options, hyp, pred_only)
n_points = size(xs,1);
% Get best set of hyperparameters
n_iter = 350;
hyp2 = [];
if ~pred_only
    evalc('[hyp2, ~, ~] = minimize(hyp, @gp, -n_iter, options.inffunc, options.meanfunc, options.covfunc, options.likfunc, x, y);');
    [~, ~, ~, ~, ~, post] = gp(hyp2, options.inffunc, options.meanfunc, options.covfunc, options.likfunc, x, y, x);
else
    hyp2 = hyp;
end

% Run the GP
if n_points > 0
    [mu, s2, fmu, fs2, lp, post] = gp(hyp2, options.inffunc, options.meanfunc, options.covfunc, options.likfunc, x, y, xs);
else
    % If just for training, return empty prediction info. The hyp2
    % structure already has the optimized values of the hyper-paramters.
    mu = []; s2 = []; fmu = []; fs2 = []; lp = [];
end
end



