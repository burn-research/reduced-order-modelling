function K = kernel(arg1, arg2, arg3, arg4)
% This function computes covariance kernels.

if nargin == 3
    K = kernel_1(arg1, arg2, arg3);
elseif nargin == 4
    K = kernel_2(arg1, arg2, arg3, arg4);
end

end

function K = kernel_1(X, type, para)
%   X: data matrix, each row is one observation, each column is one feature
%   type: type of kernel, can be 'simple', 'poly', or 'gaussian'
%   para: parameter for computing the 'poly' and 'gaussian' kernel,
%       for 'simple' it will be ignored
%   K: kernel matrix
if strcmp(type, 'simple')
    K = X*X';
elseif strcmp(type, 'poly')
    K = X*X' + 1;
    K = K.^para;
elseif strcmp(type, 'gaussian')
    K = distanceMatrix(X).^2;
    K = exp(-K./(2*para.^2)); % para is the sigma
end
end

function K = kernel_2(Y, X, type, para)
%   Y: new data matrix
%   X: training data matrix, each row is one observation, each column is one feature
%   type: type of kernel, can be 'simple', 'poly', or 'gaussian'
%   para: parameter for computing the 'poly' and 'gaussian' kernel,
%       for 'simple' it will be ignored
%   K: kernel matrix
N = size(X,1);
if strcmp(type,'simple')
    K = Y*X';
elseif strcmp(type,'poly')
    K = Y*X' + 1;
    K = K.^para;
elseif strcmp(type,'gaussian')
    K = distanceMatrix([X;Y]);
    K = K(N+1:end,1:N);
    K = K.^2;
    K = exp(-K./(2*para.^2)); % para is the sigma
end
end
