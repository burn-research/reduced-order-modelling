function [ys] = pce_predict(pol, N, c, xs, y, q, verbose, varargin)
%% Description

%% Input
if ~exist('xs', 'var') || isempty(xs)
    ys = [];
    return
end
if ~exist('y', 'var') || isempty(y)
    y = [];
end
if ~exist('verbose', 'var') || isempty(verbose)
    verbose = false;
end
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
% Predict
[Hs, ~] = pce_basis(xs, N, pol, q);
% Prediction
ys = Hs * c;
% Post-processing
if ~isempty(y)
    [y0, mu] = center(y, 1);
    [~, sig] = scale(y, y0, 2);
    ys = unscale(ys, sig);
    ys = uncenter(ys, mu);
end
end

