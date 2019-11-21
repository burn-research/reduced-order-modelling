function [Q, xp, xpnames] = getStateVars(Y, xq, vars, p, varargin)
%% Description
%{
INPUT
    Y: data matrix, rows are variables, cols are observations
    xq: spatial mesh
    vars: cell array of strings (names of the physical variables)
    p: sampled parameter

OUTPUT
    state_vars: matrix of data, rows physical variables, cols are samples
        (in the space of x and the parameter)
    xp: map returned by the pcaData function
    xpnames: cell array of the names of the parameters explored
%}

%% Default values
par_string = 'Parameter';

%% Input Analysis
% Necessary inputs
if isempty(vars)
    vars = size(Y,1);
elseif ischar(vars)
    vars = {vars};
end
if ~iscell(vars)
    % If it is not a cell, it should be a scalar (number of variables)
    vars = cell(vars, 1); % Create a cell array, whose length is the number of variables
else
    temp = vars{1};
    if ~isa(temp, 'char')
        error('Elements of the cell array VARS should be strings.');
    end
end

% Optional inputs
nin = length(varargin);
if nin > 0
    % User-supplied string (name of the parameter)
    par_string = varargin{1};
    
    % Check that the input is valid
    if ~isa(par_string,'char')
        par_string = 'Parameter';
        fprintf('User-supplied input not valid. Name of the parameter set to `Parameter`');
    end
end

%% Main
    % Useful variables
    nVars = length(vars); 
    nSamples = length(Y(1,:)); 
    n = size(xq,1);
    
    % Build state_vars
    Q = zeros(nSamples*n, nVars); % Allocate memory
    for i = 1 : nVars
        temp = [];
        for j = 1 : nSamples
            temp = [temp; Y(1+(i-1)*n:i*n,j)];
        end
        Q(:,i) = temp; % (n x nSamples) X (nVars) Data matrix
    end
    Q = Q'; % (nVars) X (n x nSamples) Data matrix
    
    % Optional output
    if nargout > 1
        xp = buildMap(xq, p);
        xpnames = {}; 
        xpnames{1} = 'x'; 
        xpnames{2} = par_string; % names of 'x' and 'p'
    end
    
end

