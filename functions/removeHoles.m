function [y] = removeHoles(x, t, val)
%% Description
% This function removes the missing number in a series of integers.
% Example: [1 2 4 5] becomes [1 2 3 4].
% Input:
% - x: (N x 1) array of integers
% - y: (N x 1) array of integers
%

%% Input
if ~exist('t', 'var') || isempty(t)
    t = false;
end
if exist('val', 'var') && ~isempty(val)
    t = true;
else
    val = 1;
end

%% Main
u = unique(x);
n_ints = length(u); % Number of unique elements in x
% Loop over values in u
for ii = n_ints - 1 : -1 : 1
    % u(ii) must equal u(ii+1) - 1
    if u(ii) + 1 ~= u(ii+1)
        m = u(ii+1) - u(ii); % Hole's size
        I = (x > u(ii)); % Values in x who are greater than u(ii)
        x(I) = (x(I) - m) + 1; % Get rid of the hole (but add 1)
    end
end
% min(x) must be 1
if t
    m = (min(x) - val);
    x = x - m;
end

%% Output
y = x;

end

