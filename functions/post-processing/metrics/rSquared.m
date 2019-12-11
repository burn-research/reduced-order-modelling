function [r2] = rSquared(Y, F, scal_crit)
%% Description
% Computes the R-squared values.
% INPUT
% - Y: (samples x variables) matrix of original data.
% - F: (samples x variables) matrix of predicted data.
% - scal_crit: scaling criterion (0 = no scaling, only centering)

%% Input
a_tol = 1e-16;
% Check they have the same size
[m, n] = size(Y);
if m ~= size(F,1) || n ~= size(F,2)
    error('Input matrices must have the same size.');
end
% Center-scale data
cent_crit = 1;
if nargin < 3
    scal_crit = 0; % Subtract the mean only
end
[Y0, ave] = center(Y, cent_crit); % Y_ave is a vector
[Y0, gamma] = scale(Y0, Y, scal_crit); % Y_gamma is a vector
[F0, ave] = center(F, cent_crit); % Y_ave is a vector
[F0, gamma] = scale(F0, F, scal_crit); % Y_gamma is a vector

%% Main
SStot = sum(Y0.^2, 1);
SSres = sum((Y0 - F0).^2, 1);

%% Output
r2 = 1 - SSres./(SStot + a_tol);

end

