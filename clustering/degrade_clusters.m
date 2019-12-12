function [new_idx] = degrade_clusters(idx)
% This function degrades cluster numeration to consecutive integers. For example, if the initial `idx` is:
%
% idx = [1,1,1,2,2,4,4,4,2]
%
% the `new_idx` returned by the `degrade_clusters()` function will be:
%
% idx = [1,1,1,2,2,3,3,3,2]
%
% so that the cluster numeration is composed of consecutive integers.
%
% Input:
% ------------
% - idx
%           a vector specifying division to clusters.
%
% Output:
% ------------
% - new_idx
%           a vector specifying division to clusters composed only of consecutive integers.

%% degrade_clusters()
% List of unique elements in idx:
u = unique(idx);

% Number of unique elements in idx:
n_ints = length(u);

% Degrade clusters:
for ii = n_ints - 1 : -1 : 1

    if u(ii) + 1 ~= u(ii+1)

        m = u(ii+1) - u(ii);
        I = (idx > u(ii));
        idx(I) = (idx(I) - m) + 1;

    end

end

% Subtract the smallest integer to make the first cluster have number 1:
m = (min(idx) - 1);
new_idx = idx - m;

end
