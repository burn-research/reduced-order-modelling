function [Yfull, x_sorted, varargout] = sortKriged(x, xk, varargin)
%% Description
%{
Given a dataset Y of observations (cols are obs) for the training points
contained in x (one column is a point), thanks to Kriging we can have
predictions for the prediction points contained in xk (one columns is a
prediction point).
This function returns a single matrix of observations Yfull, with a columns
being an observation, but the observations are sorted w.r.t. x and xk.
The sorted points are also returned in x_sorted.


INPUT
x: training points, one row is a point
xk: prediction points, one row is a point
varargin{1}: training observations, one column is an observation
varargin{2}: predicted observations, one column is an observation
varargin{3}: sorting according to the 1st, 2nd (or n-th) x dimension.
             Default is the n-th (last) dimension.


OUTPUT
Yfull: matrix of sorted observations and predictions, one column is an
observation, rows are variables.
x_sorted: sorted points.

%}



%% Input

% Dimensions
[nSamples, x_dim] = size(x);
[nkSamples, xk_dim] = size(xk);


% Varargin
nin = length(varargin);

Y = [];
Yk = [];
if nin > 0
    Y = varargin{1};
    if nin > 1
        Yk = varargin{2};
    end
end

if nin > 2
    sortingdim = varargin{3};
    if sortingdim > x_dim 
        sortingdim = x_dim; 
    end
else
    sortingdim = x_dim; % sort w.r.t. to x(:, sortingdim)
end


% Consistency on Training set
if ~isempty(Y)
    if nSamples ~= size(Y,2);
        error('Number of training points not equal to number of observations.');
    end
end

% Consistency on Prediction set
if ~isempty(Yk)
    if nkSamples ~= size(Yk,2);
        error('Number of prediction points not equal to number of observations.');
    end
end



%% Main

% Sorting order 
x_temp = x; 
if ~isempty(xk)
    x_temp = [x_temp; xk];
end
[x_temp, I] = sortrows(x_temp);



%% Output

% Sorted x
x_sorted = x_temp;

% Sorted Y
if ~isempty(Y) && ~isempty(Yk)
    Y_temp = [Y, Yk];
    Yfull = Y_temp(:,I);
elseif ~isempty(Y) && isempty(Yk)
    [~, I] = sortrows(x);
    Yfull = Y(:,I);
else
    Yfull = [];
end



end

