function pout = choosePoints(p_data, p_pred, ind)
%% Description
% Given p_data: clusters of the Local-PCA (cell array)
%       p_pred: all the unexplored points of the par-space
%       ind: # of the cluster for which we want the unexplored points
% this functions chooses which of the p2 points have to be predicted or
% explored with the p1 points corresponding to a certain cluster.
% Output: cell array


%% Main

% Number of points
nPoints = size(p_pred, 1); 

% Number of clusters
n_clust = length(p_data);

% Initialize distances
distance = zeros(n_clust,1);

% Initialize output
pout = [];


for i = 1 : nPoints
    for j = 1 : n_clust
        % Load cluster
        clust = p_data{j};
        
        % Load prediction point
        p = p_pred(i);     
        
        % Differences vector
        v = abs( clust - repmat(p, size(clust, 1), 1)); 
        
        % Distance of p(i) from clust(j)
        distance(j) = min(v); 
    end
    
    [~, j] = min(distance);
    
    if ind == j
        pout = [pout; p];
    end
end





end


