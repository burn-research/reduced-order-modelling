function [H, c, rho, k] = trainPCEKriging(X, y, npce, pol, q, cf, guess, is_parallel, verbose, varargin)
%% Description
% Trains a Kriging model that has a PCE as trend function.
%
% INPUT
% X  (obs x var): training points
% Y  (obs x var): training values
% npce: integer, number of polynomials
% cf (function handle): Kriging kernel
%
% NOTES: there is a problem with handling GUESS: should depend on x_dim,
% not y_dim

%% Input
x_dim = size(X, 2);
y_dim = size(y, 2);
% Number and degrees of polynomials
if min(size(npce)) ~= 1
    error('npce must be either a one-dim array of integers (1 x x_dim) or a scalar.');
elseif length(npce) == 1
    npce = npce * ones(1, x_dim);
elseif length(npce) ~= x_dim
    error('npce must be either a one-dim array of integers (1 x x_dim) or a scalar.');
end
% Check q (norm for hyperbolic index set)
if ~exist('q', 'var') || isempty(q)
    q = 1;
end
% Kriging kernel
if ~exist('cf', 'var') || isempty(cf)
    cf = @corrgauss; % Default 
end
% Starting guess for the evaluation of the hyperparameters for Kriging
if ~exist('guess', 'var') || isempty(guess)
    guess = [];
end
% Parallel
if ~exist('is_parallel', 'var') || isempty(is_parallel)
    is_parallel = false;
elseif ~islogical(is_parallel)
    is_parallel = (is_parallel > 0);
end
% Verbose
if ~exist('verbose', 'var') || isempty(verbose)
    verbose = false;
end
% Polynomials
pols = {'hermite', 'legendre'};
if ~exist('pol', 'var') || isempty(pol)
    pol = pols{1};
end
if ~ischar(pol)
    try
        pol = pols{round(abs(pol))};
    catch
        pol = pols{1};
    end
end

%% Main
H = cell(y_dim, 1);
c = cell(y_dim, 1);
rho = cell(y_dim, 1);
k = cell(y_dim, 1);
% Handle CF
if ~iscell(cf)
    tmp = cf;
    cf = cell(y_dim, 1);
    for ii = 1 : y_dim
        cf{ii} = tmp;
    end
end
% Handle GUESS
if length(guess) == 1 && ~isempty(guess)
    guess = zeros(y_dim, 1) + guess;
end
% Go
for ii = 1 : y_dim
    if verbose
        fprintf('Working on variable %i. \n', ii);
    end
    if isempty(guess)
        [H{ii}, c{ii}, rho{ii}, k{ii}] = local_fun(X, y(:,ii), npce, pol, q, cf{ii}, guess, is_parallel, verbose);
    else
        [H{ii}, c{ii}, rho{ii}, k{ii}] = local_fun(X, y(:,ii), npce, pol, q, cf{ii}, guess(ii), is_parallel, verbose);
    end
end
end
function [H, c, rho, k] = local_fun(X, y, npce, pol, q, cf, guess, is_parallel, verbose, varargin)
% PCE
if verbose
    fprintf('Training PCE... \n'); 
    tic;
end
[H, c, rho] = pce(X, y, npce, pol, q, verbose); 
ys = pce_predict(pol, npce, c, X, y, q, verbose); % pce_predict(H, c, X, y);
if verbose
    el_t = toc;
    fprintf('Time for PCE: %f. \n', el_t);
end
% Kriging
if verbose
    fprintf('Training Kriging... \n');
    tic;
end
tf = '';
[~, ~, k] = performKriging(X, (y - ys)', [], tf, cf, guess, is_parallel);
k = k{1};
if verbose
    el_t = toc;
    fprintf('Time for Kriging: %f. \n', el_t);
end
end
% In case of CHEBFUN
% if isa(x, 'chebfun')
%     xs = [-3:.1:3]';
%     Hs = H(xs, 1:end);
%     A = (Hs' * Hs) \ Hs';
%     ys = log(xs.^2 + 1) + 2*xs;
%     c = A * ys;
% end



