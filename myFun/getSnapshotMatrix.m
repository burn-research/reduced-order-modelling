function Y = getSnapshotMatrix(Q, xq, m, varargin)
%% Description
% INPUT
% Q: (n_observations, n_variables) data matrix
% xq: number of grid-points or grid-points
% m: number of input parameter values
%
% OUTPUT
% Y: snapshot matrix, one column of Y is one snapshot

%% Inputs
% Number of samples
if ~exist('m', 'var') || isempty(m)
    m = 1;
end
% Number of grid-points
if size(xq, 1) == 1
    nx = xq; 
else
    nx = size(xq, 1); 
end

%% Main
dim = size(Q, 2); % Number of variables
Y = zeros(dim * nx, m); % Initialize matrix
% Get snapshot matrix
for jj = 1 : m
    % Get one field, for all variables
    j2 = jj * nx;
    j1 = j2 - nx + 1;
    for ii = 1 : dim
        % Field position in snapshot matrix
        t2 = ii * nx;
        t1 = t2 - nx + 1;
        % Extract one field from Q, store it in Y
        Y(t1:t2, jj) = Q(j1:j2, ii); % Snapshot matrix
    end          
end 
end




