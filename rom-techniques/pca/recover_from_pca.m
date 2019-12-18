function [X_app_pca] = recover_from_pca(eigenvectors, scores, q, centerings, scalings)
% This function is used to recover the approximated data from Principal Component Analysis.
%
% Input:
% ------------
% - eigenvectors
%           a cell array of eigenvectors (Principal Components) in each cluster.
%
% - scores
%           a cell array of PC-scores in each cluster. It corresponds to
%           computing the equation Z = X*A' in each cluster.
%
% - q
%           the number of Principal Components (PCs) to retain in the approximation.
%           It is equivalent to specifying the rank of the approximation.
%
% - centerings
%           a vector of centerings that were applied to center the data set.
%
% - scalings
%           a vector of scalings that were applied to scale the data set.
%
% Output:
% ------------
% - X_app_lpca
%           an approximated data set.

%% recover_from_pca()
% Get dimensions:
[n_weights, n_pcs] = size(eigenvectors);

% Checks:
if ~exist('centerings', 'var') || isempty(centerings)
    centerings = zeros(1,n_weights);
end
if ~exist('scalings', 'var') || isempty(scalings)
    scalings = ones(1,n_weights);
end
if ~exist('q', 'var') || isempty(q)
    q = n_pcs;
else
  if q > n_pcs
    error('The rank requested is higher than the number of Principal Components.')
  end
end

% Recover a low-rank approximation:
X_app_pca = scores(:,1:q) * eigenvectors(:,1:q)';

% Unscale and uncenter data set:
X_app_pca = unscale(X_app_pca, scalings);
X_app_pca = uncenter(X_app_pca, centerings);

end
