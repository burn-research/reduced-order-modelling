function [B, h, R, T, P, Q, U, V] = PLS_GA(X, Y, A, varargin)
%% Description:
% INPUT
%   n X p, X
%   n X m, Y
%   number of factors, A
%

%% Input analysis
% Size check
[n, p] = size(X);
[temp, m] = size(Y);
if n ~= temp
    error('Number of rows of X must equal number of rows of Y.');
end

%% Main
% Center scale data
[X0, X_mean, X_sigma] = center_scale_data(X, 1);
[Y0, Y_mean, Y_sigma] = center_scale_data(Y, 1);
% Covariance matrix
S = X0' * Y0;
% Initialize R, T, P, Q, U and V
R = [];
T = [];
P = [];
Q = [];
U = [];
V = [];
% Start loop
for a = 1 : A
    % Dominant eigenvector of S'*S
    q = pca(S, 'Algorithm','svd', 'Centered',false);
    q = q(:,1); % [check this]
    r = S*q;
    t = X0*r;
    t = t - mean(t);
    t = t / sqrt(t'*t);
    r = r / sqrt(t'*t);
    p = X0'*t;
    q = Y0'*t;
    u = Y0*q;
    v = p;
    if a > 1
        v = v - V*(V'*p);
        u = u - T*(T'*u);
    end
    v = v / sqrt(v'*v);
    S = S - v*(v'*S);
    % Store r, t, p, q, u and v 
    R = [R, r];
    T = [T, t];
    P = [P, p];
    Q = [Q, q];
    U = [U, u];
    V = [V, v];
end

%% Output
B = R*Q';
h = diag(T*T') + 1/n;
varX = diag(P'*P)/(n-1);
varY = diag(Q'*Q)/(n-1);

end



