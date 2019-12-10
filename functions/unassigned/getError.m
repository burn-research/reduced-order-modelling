function [ER0, ER1, ER2] = getError(Y_ref, Y_app, dim, sc, varargin)
%% Description
% Input matrices Y_ref and Y_app are to be interpreted as (variables x observations).
%

%% Input
if nargin < 3
    dim = 0;
end
if nargin < 4
    sc = 0;
end

%% Main
a_tol = 1e-9;
% Get scaling factors
% [~, ~, sf] = center_scale_data(Y_ref, sc); % 1: std; 2: range
% sf = (max(Y_ref, [], 2) - min(Y_ref, [], 2));
[X0, ~, ~] = center(Y_ref', 1);
[~, sf] = scale(X0, Y_ref', sc); sf = sf(:);
% Evaluate difference matrix
ER2 = abs(Y_ref - Y_app);
% Scale the difference matrix
ER2 = ER2 .* repmat(sf, 1, size(Y_ref,2)); % Multiply
ER2 = ER2 ./ (repmat(sf.^2, 1, size(Y_ref,2)) + a_tol); % Divide
% Get a vector
ER1 = mean(ER2, 2);
% Get a scalar
ER0 = mean(ER1);

%% Output

end

