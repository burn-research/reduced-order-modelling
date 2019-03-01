function [varargout] = generateHyperparameters0(samples, values, arg1, arg2, verbose, model)
%% Description
% This script runs a lhsdesign in the length-scale space. For each point
% there, the length-scales are optimized and the corresponding surrogate
% model is judged using the LOO-analysis.
% LOO: in turn, one training observation will be left out, the
% optimized length-scales won't change, the left-out observation will be
% used as test and an error will be evaluated. The errors will be collected
% and will provide an estimate of the SM quality.
% INPUT:
% - samples, n_obs x n_var
% - values, n_var x n_obs
%

%% Input
if nargin < 5
    verbose = false;
end
if nargin < 6
    model = 'Kriging';
end
if strcmp(model, 'Kriging') || strcmp(model, 'kriging') || strcmp(model, 'k')
    is_kriging = true;
    tf = arg1;
    cf = arg2;
    dim = size(samples, 2);
else % if strcmp(model, 'GPR') || strcmp(model, 'gpr') || strcmp(model, 'Gpr')
    is_kriging = false;
    gpr_opts = arg1;
    gpr_hyp = arg2;
    dim = getHypNumber(gpr_opts, size(samples, 2));
end

%% LHS in the length-scales space
n_ls = 25; % Number of points
if is_kriging
    n_ls = 50; % Kriging is faster than GPR
end
up_ls = 100 * ones(1, dim);
lo_ls = -100 * ones(1, dim);
LS = lhsdesign_GA(n_ls, dim, up_ls, lo_ls);

%% LOO - Kriging
if is_kriging
    if verbose
        [LS_k, loo_error] = loo_kriging(samples, values, LS, tf, cf);
    else
        evalc('[LS_k, loo_error] = loo_kriging(samples, values, LS, tf, cf);');
    end
    % Get the best LS
    y_dim = size(values,1);
    ls = cell(y_dim,1);
    for ii = 1 : y_dim
        [~, minloc] = min(loo_error(:,ii));
        ls{ii} = LS_k{ii}(minloc,:);
    end
    % Output
    out_count = 0;
    if nargout > out_count
        out_count = out_count + 1;
        varargout{out_count} = ls;
    end
    return
end

%% LOO - GPR
if ~is_kriging
    if verbose
        [LS_k, loo_error] = loo_gpr(samples, values, LS, gpr_opts, gpr_hyp);
    else
        evalc('[LS_k, loo_error] = loo_gpr(samples, values, LS, gpr_opts, gpr_hyp);');
    end
    % Get the best LS
    y_dim = size(values,1);
    ls = cell(y_dim,1);
    beta = cell(y_dim,1);
    for ii = 1 : y_dim
        [~, minloc] = min(loo_error(:,ii));
        ls{ii} = LS_k{ii,1}(minloc,:);
        beta{ii} = LS_k{ii,2}(minloc,:);
    end
    % Output
    out_count = 0;
    if nargout > out_count
        out_count = out_count + 1;
        varargout{out_count} = ls;
    end
    if nargout > out_count
        out_count = out_count + 1;
        varargout{out_count} = beta;
    end
    return
end

end

function [LS_k, loo_error] = loo_kriging(samples, work_values, LS, tf, cf, varargin)
% Initialize some variables
[n_ls, x_dim] = size(LS);
y_dim = size(work_values,1); % Number of scalar targets
LS_k = cell(y_dim, 1); % To store all the matrices of optimized length-scales
loo_error = zeros(n_ls, y_dim); % To store the loo-errors
% Run LOO for each scalar target
for ik = 1 : y_dim
    fprintf('Scalar target %d out of %d\n', ik, y_dim);
    LS_opt = zeros(n_ls, x_dim); % To store the corresponding optimized values
    % Create a model for each set of candidate length-scales and evaluates
    % the LOO error (fast and slow)
    values = work_values(ik,:)'; % We need to work one target at time
    parfor ii = 1 : n_ls
        % Try to create a Kriging model. In case there's an error (probably
        % there was a problem in solving the optimization) just pretend the
        % model sucks.
        % We will find the matrix of optimized length-scales
        try
            k = get_k_model(samples, values, tf, cf, LS(ii,:));
            LS_opt(ii,:) = k.getHyperparameters();
            loo_error(ii,ik) = k.cvpe();
            % Test: try to use the gotten h-parameters as starting guess: it
            % should converge on the same starting guess. If not, there is
            % a problem!
            k_test = get_k_model(samples, values, tf, cf, k.getHyperparameters());
            hp_test = abs(k.getHyperparameters() - k_test.getHyperparameters());
            hp_test = mean(hp_test);
            if hp_test > .2
                fprintf('Test failed. Difference: %d.\n', hp_test);
            end
        catch ME
            LS_opt(ii,:) = LS(ii,:);
            loo_error(ii,ik) = Inf;
        end
    end
    LS_k{ik} = LS_opt;
end
end

function [k, m, s] = get_k_model(samples, values, trendFun, corrFun, theta0)
values = values';
evalc('[m, s, k] = performKriging(samples, values, [], trendFun, corrFun, theta0);');
end

function [LS_k, loo_error] = loo_gpr(x, values, LS, gpr_opts, gpr_hyp, varargin)
% Initialize some variables
[n_ls, dim] = size(LS);
beta_dim = length(gpr_hyp.mean);
y_dim = size(values, 1); % Number of scalar targets
LS_k = cell(y_dim, 2); % To store all the matrices of optimized length-scales
loo_error = zeros(n_ls, y_dim) + Inf; % To store the loo-errors
% Generate structures to be passed to runGpml()
hyp = cell(n_ls, 1);
m = gpr_hyp.mean;
lik = gpr_hyp.lik;
parfor ii = 1 : n_ls
    hyp{ii}.mean = m;
    hyp{ii}.cov = LS(ii,:);
    hyp{ii}.lik = lik;
end
% Run LOO for each scalar target
for ik = 1 : y_dim
    fprintf('Scalar target %d out of %d\n', ik, y_dim);
    LS_opt = zeros(n_ls, dim); % To store the corresponding optimized values
    beta_opt = zeros(n_ls, beta_dim);
    % Create a model for each set of candidate length-scales and evaluates
    % the LOO error (fast and slow)
    y = values(ik,:)'; % We need to work one target at time
    parfor ii = 1 : n_ls
        % Try to create a Kriging model. In case there's an error (probably
        % there was a problem in solving the optimization) just pretend the
        % model sucks.
        % We will find the matrix of optimized length-scales
        try
            [~, ~, hyp2] = runGpml(x, y, [], gpr_opts, hyp{ii});
            if ~isempty(hyp2.mean)
                beta_opt(ii,:) = hyp2.mean;
            end
            if ~isempty(hyp2.cov)
                LS_opt(ii,:) = hyp2.cov;
            end
            loo_error(ii,ik) = cvpe_gpml(x, y, gpr_opts, hyp2);
            % Test: try to use the gotten h-parameters as starting guess: it
            % should converge on the same starting guess. If not, there is
            % a problem!
%             k_test = get_k_model(x, y, tf, cf, k.getHyperparameters());
%             hp_test = abs(k.getHyperparameters() - k_test.getHyperparameters());
%             hp_test = mean(hp_test);
%             if hp_test > .2
%                 fprintf('Test failed. Difference: %d.\n', hp_test);
%             end
        catch ME
            if ~isempty(m)
                beta_opt(ii,:) = m;
            end
            LS_opt(ii,:) = LS(ii,:);
            loo_error(ii,ik) = Inf;
        end
    end
    LS_k{ik,1} = LS_opt;
    LS_k{ik,2} = beta_opt;
end
end

function [xval] = cvpe_gpml(x, y, options, hyp2)
[n, ~] = size(x);
xval = zeros(n,1) + Inf;
% LOO
for ii = 1 : n
    I = true(n,1);
    I(ii) = false;
    x_loo = x(I,:);
    y_loo = y(I,:);
    xs = x(~I,:);
    [mu] = gp(hyp2, options.inffunc, options.meanfunc, options.covfunc, ...
        options.likfunc, x_loo, y_loo, xs);
    xval(ii,:) = (mu - y(~I,:)).^2;
end
xval = sum(xval);
end

function dim = getHypNumber(gpr_opts, D)
try
    str = gpr_opts.covfunc();
    dim = eval(str);
catch ME
    disp(ME)
    for ii = 1 : length(ME.stack)
        disp(ME.stack(ii));
    end
    fprintf('Could not execute covfunc(). Number of hyper-parameters set to D.');
    dim = D;
end
end





