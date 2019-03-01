function [Ys, mse] = predictKriging(X, Y, Xs, k)
%% Description
% OBSERVATIONS x VARIABLES

%% Main
y_dim = size(Y, 2);
n_samples = size(Xs,1);
Ys = zeros(n_samples, y_dim);
mse = zeros(n_samples, y_dim);
for ii = 1 : y_dim
    k{ii} = k{ii}.fit(X, Y(:,ii));
    [Ys(:,ii), mse(:,ii)] = k{ii}.predict(Xs);
end
end

