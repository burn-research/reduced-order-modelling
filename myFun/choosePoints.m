function pout = choosePoints(p_data, p_pred, ind)
% Given p_data: clusters of the Local-PCA (cell array)
%       p_pred: all the unexplored points of the par-space
%       ind: # of the cluster for which we want the unexplored points
% this functions chooses which of the p2 points have to be predicted or
% explored with the p1 points corresponding to a certain cluster.
% Output: cell array

% Number of points = rows of p_pred:
nPoints=size(p_pred,1); n_clust=length(p_data); distance=zeros(n_clust,1);
pout=[];
for i=1:nPoints
    for j=1:n_clust
        clust=p_data{j};        p=p_pred(i);
        v=abs(clust-repmat(p,size(clust,1),1)); 
        distance(j)=min(v); %distance of p(i) from clust(j)
    end
    [~,j]=min(distance);
    if (ind==j)
        pout=[pout;p];
    end
end

end