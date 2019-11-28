function [r2, nrmse, mean_ab_e, rmse] = r_square(y, f, varargin)
%% Description
% Compute coefficient of determination of data fit model and RMSE
%
% [r2 nrmse] = rsquare(y,f)
% [r2 nrmse] = rsquare(y,f,c)
%
% RSQUARE computes the coefficient of determination (R-square) value from
% actual data Y and model data F. The code uses a general version of 
% R-square, based on comparing the variability of the estimation errors 
% with the variability of the original values. RSQUARE also outputs the
% root mean squared error (RMSE) for the user's convenience.
%
% Note: RSQUARE ignores comparisons involving NaN values.
% 
% INPUTS
%   Y       : Actual data
%   F       : Model fit
%
% OPTION
%   C       : Constant term in model
%             R-square may be a questionable measure of fit when no
%             constant term is included in the model.
%   [DEFAULT] TRUE : Use traditional R-square computation
%            FALSE : Uses alternate R-square computation for model
%                    without constant term [R2 = 1 - NORM(Y-F)/NORM(Y)]
%
% OUTPUT 
%   R2      : Coefficient of determination
%   NRMSE   : Normalized Root mean squared error
%   MAE     : Mean Absolute Error
%   RMSE    : Root mean squared error
%

%% Input
if isempty(varargin) 
    c = true; 
elseif length(varargin) > 1 
    error('Too many input arguments');
elseif ~islogical(varargin{1}) 
    error('C must be logical (TRUE||FALSE)');
else
    c = varargin{1};
end
% Compare inputs
if ~all( size(y) == size(f) ) 
    error('Y and F must be the same size'); 
end

ydim = size(y,2);
if ydim > 1
    r2 = zeros(ydim, 1);
    nrmse = zeros(ydim, 1);
    mean_ab_e = zeros(ydim, 1);
    rmse = zeros(ydim, 1);
    for ii = 1 : ydim
        [r2(ii), nrmse(ii), mean_ab_e(ii), rmse(ii)] = r_square(y(:,ii), f(:,ii), c);
    end
    return
end

%% Main
% Check for NaN
tmp = ~or(isnan(y),isnan(f));
y = y(tmp);
f = f(tmp);
if c 
    r2 = max(0, 1 - sum((y(:) - f(:)).^2) / sum((y(:) - mean(y(:))).^2));
else
    r2 = 1 - sum((y(:) - f(:)).^2) / sum(y(:).^2);
    if (r2 < 0)
    % http://web.maths.unsw.edu.au/~adelle/Garvan/Assays/GoodnessOfFit.html
        warning('Consider adding a constant term to your model') %#ok<WNTAG>
        r2 = 0;
    end
end
nrmse = sqrt(mean((y(:) - f(:)).^2)) / sqrt(mean(y.^2));
mean_ab_e = mean(abs(y(:) - f(:))) / mean(abs(y));
rmse = sqrt(mean((y(:) - f(:)).^2));
end


