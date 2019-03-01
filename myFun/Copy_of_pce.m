function [H, c, rho] = pce(X, y, N, pol, varargin)
%% Description
% Trains a SM based on Polynomial Chaose Expansion (PCE).
%
% INPUT
% X (M x x_dim) are the training locations.
% y (M x 1) are the training values.
% N (1 x x_dim) are the degrees of the polymonials for each input variable.
%

%% Input
% Check input type
if ~isnumeric(X) || ~isnumeric(y)
    error('Provide numeric data.');
end
x_dim = size(X, 2); % Input dimension
% Check validity for N
if min(size(N)) ~= 1
    error('N must be either a one-dim array of integers (1 x x_dim) or a scalar.');
elseif length(N) == 1
    N = N * ones(1, x_dim);
elseif length(N) ~= x_dim
    error('N must be either a one-dim array of integers (1 x x_dim) or a scalar.');
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
% Pre-process data
y0 = center(y, 1);
y = scale(y, y0, 3); clear y0
% Start
lb = []; ub = [];
if x_dim > 1
    lb = min(min(X, [], 1));
    ub = max(max(X, [], 1));
end
H = cell(x_dim,1);
c = cell(x_dim,1);
rho = cell(x_dim,1);
for ii = 1 : x_dim
    [H{ii}, c{ii}, rho{ii}] = pce_1(X(:,ii), y, N(ii), lb, ub, pol);
end
% Output
if x_dim == 1
    H = H{1};
    c = c{1};
    rho = rho{1};
else
    % Product of univariate polynomials 
    H_tmp = [];
    for i = 1 : x_dim - 1
        for j = i + 1 : x_dim
            for k = 1 : N(i) + 1
                for l = 1 : N(j) + 1
                    h = H{i}(X(:,i),k) .* H{j}(X(:,j),l);
                    H_tmp = [H_tmp, h];
                end
            end
        end
    end
    % Determination of coefficients
    A = (H_tmp' * H_tmp) \ H_tmp';
    c = A * y; % Coefficients
    mse = immse(y, H_tmp * c);
    % Lasso
    [c_lass, info] = lasso(H_tmp, y);
    % Choose c
    [~, idx] = knee_pt([mse, info.MSE]);
%     [~, idx] = min([mse, info.MSE]);
    c = [c(:), c_lass];
    c = c(:,idx);
end
end


function [H, c, rho] = pce_1(x, y, N, lb, ub, pol, varargin)
%% Description
% NOTES
% Legendre:
% P(n,x) = (2*n-1)/n * x * P(n-1,x) - (n-1)/n * P(n-2,x)
%
% Should data be centered here ?

%% Input
if ~isnumeric(x) || ~isnumeric(y)
    error('Provide numeric data.');
end
if min(size(x)) ~= 1 || min(size(y)) ~=1
    error('Provide scalar data: X and Y must be vectors.');
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
% Column vectors
x = x(:); 
y = y(:);
% Range for x
if ~exist('lb', 'var') || ~exist('ub', 'var') || isempty(lb) || isempty(ub)
    lb = min(x);
    ub = max(x);
end
% Hermite polynomials
z = chebfun(@(x) x, [lb ub]);
H = [z*0 + 1, z];
switch pol
    case pols{1}
        % Hermite
        for n = 3 : N + 1
            H(:,n) = z.*H(:,n-1) - n*H(:,n-2);
%             H(:,n) = H(:,n) / norm(H(:,n));
        end
    case pols{2}
        % Legendre
        for n = 3 : N + 1
            H(:,n) = (2*n-1)/n * z.*H(:,n-1) - (n-1)/n * H(:,n-2);
%             H(:,n) = H(:,n) / norm(H(:,n));
        end
end
% Gaussian density
rho = exp(-z.^2 / 2);  
rho = rho / sum(rho);
% Determination of coefficients
A = (H(x,:)' * H(x,:)) \ H(x,:)';
c = A * y; % Coefficients
mse = immse(y, H(x,:) * c);
% Lasso
[c_lass, info] = lasso(H(x,:), y);
% Choose c
[~, idx] = knee_pt([mse, info.MSE]);
% [~, idx] = min([mse, info.MSE]);
c = [c(:), c_lass];
c = c(:,idx);
end


