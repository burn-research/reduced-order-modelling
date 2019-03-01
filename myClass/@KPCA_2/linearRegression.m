function [predictions, mse, k] = linearRegression(samples, values, prediction_points, varargin)
%% INPUTS:
%    x - independent variables.  These should be ordered with each
%        observation in a separate row and each variable in a column.
%    y - dependent variables.  These should be ordered with each observation
%        in a separate row and each variable in a column.
%   OPTIONAL INPUTS:
%    atol    - refinement tolerance
%    maxOrd  - maximum order for basis functions (default is 3)
%    nxtrial - number of (uniform) bins in x space to use (default is 4)
%    maxiter - maximum number of refinement iterations (default is 30)
%   

%% Options
tol = 1e-5;
if ~isempty(varargin)
    maxOrd = varargin{1};
else
    maxOrd = 3;
end
nxtrial = 4;
maxiter = 50;

%% MARS
% Initialize variables
predictions = zeros( size(values,1), size(prediction_points,1) );
mse = predictions * 0;
k = cell(size(values, 1),1);

% Run MARS
parfor i = 1 : size(values, 1)
    temp_model = mars(samples, values(i,:)', tol, maxOrd, nxtrial, maxiter);
    predictions(i,:) = temp_model(prediction_points)';
    k{i} = temp_model;
end


end


