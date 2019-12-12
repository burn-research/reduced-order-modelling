function [valsout, varargout] = Krigeit(samples, values, prediction_point, varargin)
%% Ddescription
%{
INPUT
samples: one row is a sample, the dimension of the parameter is the number of cols
values: one row is an observation, cols is the dimension of the observation

OUTPUT
valsout: predictions, one row is a prediction
%}



%% Default
trendFun = 'regpoly2'; %default trend function
corrFun = @corrmatern32; %default corr function
% corrFun=@corrgauss;



%% Input

% Dimensions
[nSamples, dim] = size(samples);
d = dim;
nin = length(varargin);

if ~isempty(varargin) %check further inputs 
    trendFun = varargin{1};
    corrFun = varargin{2};
end 

if ~strcmp(trendFun(1:end-1),'regpoly') && ~strcmp(trendFun,'')
    fprintf('\tRegression function was mispecified, the default one was set.\t');
    trendFun='regpoly2'; %if the trend fun was not properply specified, choose the default one
end

% Hyperparameters guess
val = .1;
if nin > 2
    val = varargin{3};
end


%% Options
opts = Kriging.getDefaultOptions();
opts.hpOptimizer = SQPLabOptimizer( dim, 1 );
opts.regressionMaxOrder = dim;



%% Main 

% Hyperparameters 
Val = 1e10;
lb = zeros(1,d);  ub = Val * ones(1,d); % 1 x d
theta0 = val * ones(1,d);               % 1 x d
opts.hpBounds = [lb; ub];               % 2 x d hyperparameter optimization bounds


% Kriging constructor
k = Kriging( opts, theta0, trendFun, corrFun );
k = k.fit(samples, values);


% Predict
[temp, temp2] = k.predict(prediction_point); 

% temp = temp';   temp2 = temp2';

% NaN check
checkNaN = isnan(temp);   temp(checkNaN) = 0;
checkNaN = isnan(temp2);  temp2(checkNaN) = 0;


%% Output
valsout = temp;
mse = temp2;



%% Varargout
if nargout > 1
    varargout{1} = mse;
end



end




