function D = distanceMatrix(X)
%   X: data matrix, each row is one observation, each column is one feature
%   D: pair-wise distance matrix
%
N = size(X, 1);
XX = sum(X.*X, 2);
D = repmat(XX, 1, N) + repmat(XX', N, 1) - 2*(X*X');
D(D<0) = 0;
D = sqrt(D);
end