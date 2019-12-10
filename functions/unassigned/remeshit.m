function [val, mse] = remeshit(X, Y, prediction_points, varargin)
%% Description
%{
Works for 2D meshes (3rd method works with 1D mesh). 

INPUTS
    X: cell array. X{i} is the N-by-2 mesh for the i-th case.
    Y: cell array. Y{i} is the N-by-V mesh for the i-th case, V can be any.
    prediction_points: common mesh, M-by-2 matrix, with M being any.
    varargin{1}: method: Kriging() or scatteredInterpolant().

OUTPUTS
    val: interpolated values, val{i} is a M-by-V matrix.
    mse: mean squared errors, only for the Kriging() case.

%}


%% Main

% Additional inputs
nin = length(varargin);     


% Check for the method
method = 2; % Default
if nin > 0
    method = varargin{1};
end


% Initialize outputs
val = cell(size(Y));    mse = val;



% case 1 (KRIGING)
if method == 1
    for i = 1 : length(Y)
        if isa(X,'cell')
            samples = X{i};     
        else
            samples = X;
        end
        values = Y{i};

        [val{i}, mse{i}] = Krigeit(samples, values, prediction_points);
    end
    
    return;
end



% case 2 (GRIDDATA)
if method == 2
    % Separate the columns
    xq = prediction_points(:,1);
    yq = prediction_points(:,2);
    
    % Go case by case
    for i = 1 : length(Y)
        if isa(X,'cell')
            samples = X{i};     
        else
            samples = X;
        end
        x = samples(:,1);  y = samples(:,2);
        values = Y{i};
    
        % Variable by variable, for the same case
        for j = 1 : size(values, 2)
            F = scatteredInterpolant(x, y, values(:,j));
            val{i}(:,j) = F(xq, yq);
        end
    end 
    
    return;
end 




% case 3 (1D mesh)
if method == 3
    % Check everything is fine
    if size(prediction_points,1) ~= 1 && size(prediction_points,2) ~= 1 
        error('PREDICTION_POINTS should be a 1D mesh if METHOD is 3.');
    end
    
    % Separate the columns
    xq = prediction_points;
    interp_method = 'spline';
    
    % Go case by case
    for i = 1 : length(Y)
        if isa(X,'cell')
            samples = X{i};     
        else
            samples = X;
        end     
        x = samples(:,1);
        values = Y{i};
    
        % Variable by variable, for the same case
        for j = 1 : size(values, 2)
            val{i}(:,j) = interp1(x, values(:,j), xq, interp_method);
        end
    end 
    
    return;
end 


end


%       'nearest'  - nearest neighbor interpolation
%       'next'     - next neighbor interpolation
%       'previous' - previous neighbor interpolation
%       'linear'   - linear interpolation
%       'spline'   - piecewise cubic spline interpolation (SPLINE)
%       'pchip'    - shape-preserving piecewise cubic interpolation
%       'cubic'    - same as 'pchip'
%       'v5cubic'  - the cubic interpolation from MATLAB 5, which does not
%                    extrapolate and uses 'spline' if X is not equally
%                    spaced.

