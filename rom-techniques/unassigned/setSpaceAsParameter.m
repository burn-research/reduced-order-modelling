function [xp, I] = setSpaceAsParameter(xy, p)
%% Description
%{
INPUT
    p: matrix of parameters values. Rows are samples, columns are number of
       parameters.
    xy: mesh. Rows are number of grip-points, columns are dimensionality of
        the mesh.

OUTPUT
    xp: combined matrix of parameters and grid-points.
%}




%% Case where the mesh is not common for every case (must be a cell array)
if isa(xy,'cell')
    [xp, I] = cellArrayCase(xy, p);
    return;
end




%% Case where the mesh is common for every case

% Sizes
[nSamples, nParameters] = size(p);
[nGridPoints, meshDimension] = size(xy);


% Initialize the output
A = zeros(nSamples * nGridPoints, nParameters + meshDimension);


for i = 1 : nSamples
    % Repeat the parameter set
    P = repmat(p(i,:), nGridPoints, 1); 
    
    % Get the starting and final indeces
    i_in = 1 + nGridPoints * (i-1);
    i_f = i * nGridPoints;
    
    % Fill the output matrix
    A(i_in:i_f, :) = [xy, P];
end


% Sorting
sortingdim = meshDimension + nParameters; % Choose the column to sort
x_temp = A(:, sortingdim); % Get that column
[~, I] = sort(x_temp); % Get the sorted indeces


%% Output
xp = A(I,:);




end


function [xp, I] = cellArrayCase(xy, p)

% Notes: numberOfMeshes = nSamples

% Sizes
[nSamples, nParameters] = size(p);
numberOfMeshes = length(xy);


% Check
if numberOfMeshes ~= nSamples
    error('Number of samples must equal the number of meshes.');
end


% Initialize the output
A = [];

for i = 1 : nSamples
    % Sizes
    [nGridPoints, meshDimension] = size(xy{i});
    
    % Repeat the parameter set
    P = repmat(p(i,:), nGridPoints, 1); 
    
    % Fill the output matrix
    temp = [xy{i}, P];      A = [A; temp];
end


% Sorting
sortingdim = meshDimension + nParameters; % Choose the column to sort
x_temp = A(:, sortingdim); % Get that column
[~, I] = sort(x_temp); % Get the sorted indeces


%% Output
xp = A(I,:);


end







