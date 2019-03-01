% Get the data correspongind to the value of the parameters equal to xp.
% INPUTS:
%           x  = the point
%           xp = the map
%           Y  = the data (datadim X length of the map)

function value = sampleit(x,xp,Y)
[~,indx] = ismember(x,xp,'rows');

value = Y(:,indx);
end