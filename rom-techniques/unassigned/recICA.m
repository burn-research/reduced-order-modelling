function [Y_ica] = recICA(varargin)
%% Input
if nargin == 2
    W = varargin{1};
    Y0 = varargin{2};
elseif nargin == 3
    W = varargin{1};
    S = varargin{2};
    T = varargin{3};
end
%% Main
if exist('S', 'var')
    Y_ica = (S * W') / T'; % [case n x r]
else
    [Z, T] = whitenRows(Y0');
    T = T'; Z = Z';
    S = Z * W;
    Y_ica = (S * W') / T; % [case n x r]
end
end

% 
% [Z, T] = whitenRows(Y0');
% S = W * Z;
% S = (W' * Z)';

