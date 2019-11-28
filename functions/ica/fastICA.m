function [S, W, T, nrmse] = fastICA(Z, r, is_white, type, flag, TOL, MAX_ITERS)
%% Description
% Syntax:       S = fastICA(Z,r);
%               S = fastICA(Z,r,type);
%               S = fastICA(Z,r,type,flag);
%               [S, W, T, mu] = fastICA(Z,r);
%               [S, W, T, mu] = fastICA(Z,r,type);
%               [S, W, T, mu] = fastICA(Z,r,type,flag);
%               
% Inputs:       Z is an n x d matrix containing n samples of d-dimensional
%               data
%               
%               r is the number of independent components to compute
%
%               [OPTIONAL] is_white, will whiten data if TRUE
%               
%               [OPTIONAL] type = {'kurtosis','negentropy'} specifies
%               which flavor of non-Gaussianity to maximize. The default
%               value is type = 'negentropy'
%               
%               [OPTIONAL] flag determines what status updates to print
%               to the command window. The choices are
%                   
%                       flag = 0: no printing
%                       flag = 1: print iteration status
%               
% Outputs:      S is an n x r matrix containing the r independent
%               components - scaled to variance 1 - of the input samples
%               
%               W and T are the ICA transformation matrices such that
%               Zr = T \ (W' * S) + repmat(mu,1,n);
%               [Zr = S * W'] (Added by me)
%               is the r-dimensional ICA approximation of Z
%               
%               mu is the d x 1 sample mean of Z
%               
% Description:  Performs independent component analysis (ICA) on the input
%               data using the Fast ICA algorithm
%               
% Reference:    Hyv?rinen, Aapo, and Erkki Oja. "Independent component
%               analysis: algorithms and applications." Neural networks
%               13.4 (2000): 411-430
%               
% Author:       Brian Moore
%               brimoor@umich.edu               
% Date:         April 26, 2015
%               November 12, 2016
%               May 4, 2018
%
% Edit:         Gianmarco Aversano, Gianmarco.Aversano@ulb.ac.be
% Date:         October, 2018
%
%% Input
% Parse inputs
Z = Z';
if ~exist('is_white','var') || isempty(is_white)
    % Default is_white
    is_white = true;
end
if ~exist('flag','var') || isempty(flag)
    % Default display flag
    flag = true;
end
if ~exist('type','var') || isempty(type)
    % Default type
    type = 'negentropy';
end
if ~exist('TOL','var') || isempty(TOL)
    % Default TOL
    TOL = 1e-6; % Convergence criteria
end
if ~exist('MAX_ITERS','var') || isempty(MAX_ITERS)
    % Default MAX_ITERS
    MAX_ITERS = 100; % Max # iterations
end
n = size(Z, 1);
a_tol = 1e-16;
%% Set algorithm type
if strncmpi(type,'kurtosis',1)
    % Kurtosis
    USE_KURTOSIS = true;
    algoStr = 'kurtosis';
elseif strncmpi(type,'negentropy',1)
    % Negentropy
    USE_KURTOSIS = false;
    algoStr = 'negentropy';
else
    % Unsupported type
    error('Unsupported type ''%s''',type);
end
%% Main
% Center and whiten data
% [Zc, mu] = centerRows(Z);
if is_white
    [Z, T] = whitenRows(Z);
else
    T = 1;
end
% Normalize rows to unit norm
normRows = @(X) bsxfun( @rdivide, X, sqrt(sum(X.^2, 2) + 1e-16) );
% Perform Fast ICA
if flag
    % Prepare status updates
    fmt = sprintf('%%0%dd', ceil(log10(MAX_ITERS + 1)));
    str = sprintf('Iter %s: max(1 - |<w%s, w%s>|) = %%.4g\\n', fmt, fmt, fmt);
    fprintf('***** Fast ICA (%s) *****\n',algoStr);
end
W = normRows(rand(r, size(Z, 1))); % Random initial weights
k = 0;
delta = inf;
while (delta > TOL) && (k < MAX_ITERS)
    k = k + 1;
    % Update weights
    Wlast = W; % Save last weights
    Sk = W * Z;
    if USE_KURTOSIS
        % Kurtosis
        G = 4 * Sk.^3;
        Gp = 12 * Sk.^2;
    else
        % Negentropy
        G = Sk .* exp(-0.5 * Sk.^2);
        Gp = (1 - Sk.^2) .* exp(-0.5 * Sk.^2);
    end
    W = (G * Z') / n - bsxfun(@times,mean(Gp,2),W);
    W = normRows(W);
    % Decorrelate weights
    [U, S, ~] = svd(W, 'econ');
    W = U * diag(1 ./ diag(S + a_tol)) * U' * W;
    % Update convergence criteria
    delta = max(1 - abs(dot(W, Wlast, 2)));
    if flag
        fprintf(str, k, k, k - 1, delta);
    end
end
if flag
    fprintf('\n');
end
% Independent components
S = W * Z; % S = Z * W;
% Zr = T \ (W' * S); % [case r x n]
% Zr = (S * W') / T; % [case n x r]
%% Output
S = S';
W = W';
T = T';
% Reconstruct and get errors
Z = Z';
Z_ica = (S * W') / T;
[~, nrmse, ~, ~] = r_square(Z, Z_ica);
end



