function value = sampleit(x, xp, Y)
% Get the data correspongind to the value of the parameters equal to xp.
% INPUTS:
%           x  = the point
%           xp = the map
%           Y  = the data (datadim X length of the map)

[~, indx] = ismember(x, xp, 'rows');
value = Y(:, indx);

end