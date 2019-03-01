function [idx, uncentered_nz_X_k, this] = localPCA_hierarchical(X, n_eigs, cutoff, maxclust)


% Dimensions
[rows, columns] = size(X);
n_eigs_max = columns;


% Clustering options
distance = @distfun;


% Cluster data
idx = clusterdata(X, ...
    'criterion', 'distance', ...
    'cutoff', cutoff,...
    'distance', distance, ...
    'maxclust', maxclust);


% Output
uncentered_nz_X_k = []; 
this = [];

end


function d2 = distfun(XI, XJ)
% DESCRIPTION
%     T = clusterdata(X, cutoff);
%     T = clusterdata(X, Name, Value);
% 
% A distance function specified using @:
%     D = pdist(X, @distfun);
% 
% A distance function must be of the form:
%     d2 = distfun(XI, XJ);
%     
%     XI: (1 x n) single row of X
%     XJ: (m2 x n) multiple rows of X (arbitrary number of rows)
%     d2: (m2 x 1) distances, whose h-th elements is the distance between XI
%         and XJ(k,:)
% 

% Size
[rows, cols] = size(XJ);

% Apply PCA
[modes, scores, ~] = pca(XJ');  
scores = scores'; 

% Approximation order
n_eigs = cols - 2;

% Recover XI
XI_rec = modes * scores;

% Distances
d2 = .5 * norm(XI - XI_rec);

end





