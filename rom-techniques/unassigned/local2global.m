function [A, varargout] = local2global(firstCell, secondCell, varargin)
%{
DESCRIPTION
From a set of local data, rebuild the global data.

INPUT
    firstCell: each cell is a matrix of local data.
    secondCell: each cell is a set of row-indeces, used to sort the local data.

OUTPUT
    A: global matrix
%}


%% Input

% Check inputs are CELL arays
if ~isa(firstCell,'cell') || ~isa(secondCell,'cell')
    error('Input must be CELL arrays');
end

% Dimensions
numberOfClusters = length(firstCell);
if numberOfClusters ~= length(secondCell);
    error('Dimensions of inputs mismatch.');
end

% varargin
nin = length(varargin);
s = [];
if nin > 0
    s = varargin{1};
    if ~isa(s,'char')
        error('Third input must be a CHAR.');
    end
end

% Check class of each cell
if ~isa(firstCell{1},'double') 
    prompt = 'Which field to look for?\n';
    if isempty(s)
        s = input(prompt,'s');
    end
    
    for i=1:numberOfClusters
        temp = eval( ['firstCell{i}.',s] );
        firstCell{i} = temp;
    end
    
end


%% Main

% Get the number of rows and cols
temp = 0;
for i=1:numberOfClusters
    temp = temp + size(firstCell{i}, 1);
end
rows = temp; 
cols = size(firstCell{1}, 2);

% Initialize A
A_out = zeros(rows, cols);

% Fill A
for k=1:numberOfClusters
    temp = secondCell{k}; % Load cluster
    for j=1:length(temp) 
        i = temp(j); % Get rows associated with this cluster
        A_out(i,:) = firstCell{k}(j,:); % Write them on the right position inside A
    end
end


%% Output
A = A_out;


end


