function [ys, mse] = predictPCEKriging(Xs, pol, N, c, q, k, Y, verbose, varargin)
%% Description
%
% INPUT
%

%% Input
if ~exist('verbose', 'var') || isempty(verbose)
    verbose = false;
end
% Check q (norm for hyperbolic index set)
if ~exist('q', 'var') || isempty(q)
    q = 1;
end

%% Main
if iscell(k)
    y_dim = length(k);
    ys = zeros(size(Xs,1), y_dim);
    mse = zeros(size(Xs,1), y_dim);
    for ii = 1 : y_dim
        % PCE
        ypce = pce_predict(pol, N, c{ii}, Xs, Y(:,ii), q, verbose); % pce_predict(H, c, Xs, Y);
        % Kriging
        [yk, mse] = k{ii}.predict(Xs);
        % Output
        ys(:,ii) = ypce + yk;
    end
else
    % PCE
    ypce = pce_predict(pol, N, c, Xs, Y, q, verbose); % pce_predict(H, c, Xs, Y);
    % Kriging
    [yk, mse] = k.predict(Xs);
    % Output
    ys = ypce + yk;
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



