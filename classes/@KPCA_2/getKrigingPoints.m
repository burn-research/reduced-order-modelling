function getKrigingPoints(obj)
%% Description
% Given p_data: clusters of the Local-PCA (cell array)
%       p_pred: all the unexplored points of the par-space
%       ind: # of the cluster for which we want the unexplored points
% this function chooses which of the p2 points have to be predicted or
% explored with the p1 points corresponding to a certain cluster.
% Output: cell array


%% Main

% Number of candidate Kriging points
nPoints = size(obj.prediction_points, 1); 

% Initialize distances
distance = zeros(obj.number_of_clusters, 1);

for i = 1 : nPoints
    
    for j = 1 : obj.number_of_clusters        
        % Load prediction point
        p = obj.prediction_points(i, :);     
        
        % Euclidean distances of p from each point of obj.lpca{j}.xp
        cols = size(obj.local_pca{j}.training_data, 2);
        v = obj.local_pca{j}.training_points - repmat(p, cols, 1); 
        v = sum(v.^2, 2).^.5;
        
        % Distance of p from cluster j
        distance(j) = min(v); 
    end
    
    % Evaluate the cluster to which P is the closest
    [~, j] = min(distance);
    
    % Append that point
    xp_kriged = [obj.local_pca{j}.prediction_points; p];
    obj.lpca{j}.prediction_points = xp_kriged;
    
end


end


