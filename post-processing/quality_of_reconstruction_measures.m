function [r2, nrmse, mae, rmse] = quality_of_reconstruction_measures(X, F, r2_constant)
% This function computes various quality of reconstruction measures such as
% the coefficient of determination (R2), Normalized Root Mean Squared Error (NRMSE),
% Mean Absolute Error (MAE) and Root Mean Squared Error (RMSE).
%
% Input:
% ------------
% - X
%     the original data set.
%
% - F
%     a model fit to the original data set.
%
% - r2_constant
%     optional Boolean specifying whether R2 computation
%     according to constant term in the model:
%     https://web.maths.unsw.edu.au/~adelle/Garvan/Assays/GoodnessOfFit.html
%     If this value is not specified or set to false, a traditional R2 computation is performed.
%
% Output:
% ------------
% - r2
%     coefficient of determination.
%
% - nrmse
%     Normalized Root Mean Squared Error.
%
% - mae
%     Mean Absolute Error.
%
% - rmse
%     Root Mean Squared Error.

%% quality_of_reconstruction_measures()
% Checks:
if ~all( size(X) == size(F) )
    error('Y and F must be the same size');
end

if ~exist('r2_constant', 'var') || isempty(r2_constant)
    r2_constant = false;
end

ydim = size(X,2);

if ydim > 1

    % Initialize measures:
    r2 = zeros(ydim, 1);
    nrmse = zeros(ydim, 1);
    mae = zeros(ydim, 1);
    rmse = zeros(ydim, 1);

    % Run this function in recurrence on each column of X and F:
    for ii = 1:1:ydim
        [r2(ii), nrmse(ii), mae(ii), rmse(ii)] = quality_of_reconstruction_measures(X(:,ii), F(:,ii), r2_constant);
    end

    return
end

% Coefficient of determination:
if r2_constant == false
  r2 = (1 - sum((X(:) - F(:)).^2) / sum((X(:) - mean(X(:))).^2));
else
  % Computation according to constant term in the model:
  r2 = 1 - sum((X(:) - F(:)).^2) / sum(X(:).^2);
  if (r2 < 0)
      warning('Consider adding a constant term to your model.')
      r2 = 0;
  end
end

% Normalized Root Mean Squared Error:
nrmse = sqrt(mean((X(:) - F(:)).^2)) / sqrt(mean(X.^2));

% Mean Absolute Error:
mae = mean(abs(X(:) - F(:))) / mean(abs(X));

% Root Mean Squared Error:
rmse = sqrt(mean((X(:) - F(:)).^2));

end
