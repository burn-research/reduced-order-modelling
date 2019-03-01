function [X, Y, Z] = myMeshGrid(x, y, z)

% Input
if min(size(x)) ~= 1 || min(size(y)) ~= 1 || min(size(z)) ~= 1
    error('Inputs must be vectors.');
end
if length(x) ~= length(y) || length(x) ~= length(z)
    error('Inputs must have same length.');
end

% Make them column vectors
x = x(:);   
y = y(:);   
z = z(:);


% Size
l = length(x);


% Check if there are doubles
xy = [x,y];
if size(unique(xy,'rows'), 1) ~= l
    error('There are repeated points.');
end


% Default meshgrid function
[X, Y] = meshgrid(x, y);


% Z
Z = zeros(size(X));
for i = 1 : l
    for j = 1 : l
        xi = X(j,i);     yj = Y(j,i);
        
        temp = xy - repmat([xi,yj], l, 1);
        [row, col] = find(temp);
        
        thisz = z(row(1));
        
        Z(j,i) = thisz;
    end
end


end



% [X,Y] = meshgrid(1:3,10:14)
% X =
%      1     2     3
%      1     2     3
%      1     2     3
%      1     2     3
%      1     2     3
% Y =
%     10    10    10
%     11    11    11
%     12    12    12
%     13    13    13
%     14    14    14


