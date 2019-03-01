function [H, rho, A] = pce_basis(X, N, pol, q, varargin)
%% Description

%% Input
% Sizes
[L, x_dim] = size(X);
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
% Check validity for N
if min(size(N)) ~= 1
    error('N must be either a one-dim array of integers (1 x x_dim) or a scalar.');
elseif length(N) == 1
    N = N * ones(1, x_dim);
elseif length(N) ~= x_dim
    error('N must be either a one-dim array of integers (1 x x_dim) or a scalar.');
end
% Check q (norm for hyperbolic index set)
if ~exist('q', 'var') || isempty(q)
    q = 1;
end

%% Main
A = [];
if x_dim == 1
    % Univariate
%     H = legendre(N, X, 'norm')';
    H = zeros(L, N+1);
    for ii = 0 : N
        tmp = legendre(ii, X, 'norm')';
        H(:,ii+1) = tmp(:,1);
    end
    % Gaussian density
    rho = exp(-X.^2 / 2);  
    rho = rho / sum(rho);
else
    % Multivariate
    Huni = cell(x_dim, 1);
    for ii = 1 : x_dim
%         Huni{ii} = legendre(N(ii), X(:,ii), 'norm')';
        Huni{ii} = zeros(L, N(ii)+1);
        for jj = 0 : N(ii)
            tmp = legendre(jj, X(:,ii), 'norm')';
            Huni{ii}(:,jj+1) = tmp(:,1);
        end
    end
    % Get hyperbolic index set
    p = min(N);
    [A] = index_set(x_dim, p, q);
    % Build multivariate polynomials
    H = [];
    for ii = 1 : size(A, 1)
        h = 1;
        for jj = 1 : x_dim
            h = h .* Huni{jj}(:,A(ii,jj));
        end
        H = [H, h];
    end
    % Gaussian density
    tmp = exp(-X.^2 / 2); 
    rho = ones(size(X,1), 1);
    for ii = 1 : x_dim
        rho = rho .* tmp(:,ii);
    end
    rho = rho / sum(rho);
end
end


function [H, rho] = keep_code(varargin)
if x_dim == 1
    % Univariate
    X = X(:); % column-vector
    % Polynomials
    H = size(m, N);
    for k = 1 : N
        switch pol
            case pols{1}
                h = hermiteH(k, X);
            case pols{2}
                h = legendreP(k, X);
        end
        H(:,k) = h;
    end
    % Gaussian density
    rho = exp(-x.^2 / 2);  
    rho = rho / sum(rho);
else
    % Multivariate
    H = [];
    it = 0;
    for i = 1 : x_dim - 1
        for j = i + 1 : x_dim
            for k = 1 : N(i) + 1
                for l = 1 : N(j) + 1
                    switch pol
                        case pols{1}
                            h = hermiteH(k, X(:,i)) .* hermiteH(l, X(:,j));
                        case pols{2}
                            h = legendreP(k, X(:,i)) .* legendreP(l, X(:,j));
                    end
                    it = it + 1;
                    H = [H, h];
                end
            end
        end
    end
    rho = []; % returned empty
end
% more
H = [];
    p = min(N) + 1;
    for i = 1 : x_dim - 1
        for j = i + 1 : x_dim
            for k = 1 : N(i) + 1
                for l = 1 : N(j) + 1
                    
                    h = Huni{i}(:,k) .* Huni{j}(:,l);
                    H = [H, h];
                end
            end
        end
    end
end
