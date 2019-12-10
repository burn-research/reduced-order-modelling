function [H, c, rho, Ai] = pce(X, y, N, pol, q, verbose, varargin)
%% Description
% Trains a SM based on Polynomial Chaose Expansion (PCE).
%
% INPUT
% X (M x x_dim) are the training locations.
% y (M x 1) are the training values.
% N (1 x x_dim) are the degrees of the polymonials for each input variable.
%

%% Input
% Verbose
if ~exist('verbose', 'var') || isempty(verbose)
    verbose = false;
end
% Check input type
if ~isnumeric(X) || ~isnumeric(y)
    error('Provide numeric data.');
end
% Polynomials
pols = {'hermite', 'legendre'};
if ~exist('pol', 'var') || isempty(pol)
    pol = pols{2};
end
if ~ischar(pol)
    try
        pol = pols{round(abs(pol))};
    catch
        pol = pols{2};
    end
end
% Check q (norm for hyperbolic index set)
if ~exist('q', 'var') || isempty(q)
    q = 1;
end

%% Main
% Pre-process data
y0 = center(y, 1);
y = scale(y, y0, 2); clear y0
% Get H
[H, rho, Ai] = pce_basis(X, N, pol, q);
% Determination of coefficients
A = (H' * H) \ H';
c = A * y; % Coefficients
mse = immse(y, H * c);
% Lasso
[c_lass, info] = lasso(H, y, 'Alpha', 1e-3); % 'CV', 4,
% Choose c
[~, idx] = knee_pt([mse, info.MSE]);
% [~, idx] = min([mse, info.MSE]);
c = [c(:), c_lass];
c = c(:,idx);
end
