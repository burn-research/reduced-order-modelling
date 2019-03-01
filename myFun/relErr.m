function err_rel = relErr(Y, Y_full, varargin)
%{
This function computes the relative error between the matrices or vectors 
Y and Y_full (differences between columns).
In case of 2 matrices, the 3rd input is a normalization factor.
In case of 2 vectors, the 3rd input is the flag for how to treat the zeros.
%}

%% Input Check

% Centered-scaled data
% [Y, ave, sig] = zscore(Y, 0, 2);
% [Y_full, ave_full, sig_full] = zscore(Y_full, 0, 2);


% Sizes
[m, n] = size(Y);


% Check dimensions are equal
if m ~= size(Y_full,1) || n ~= size(Y_full,2)
    error('You`re trying to compare variables with different sizes.');
end


% Check if input is a vector or matrix
isVector = false;
if m == 1 || n == 1
    isVector = true;
end


% Normalization factor
norm_factor = zeros(1, n);
for i = 1 : n
    norm_factor(:,i) = norm(Y_full(:,i)) + eps; % row vector
end


% Zero Strategy
zeroStrategy = false;


% Get additional inputs
third_input = [];
if nargin > 2
    third_input = varargin{1};
    if ~isVector && ~isa(third_input,'double')
        error('Something is wrong.');
    end
    if isVector && ~isa(third_input,'logical')
        error('Something is wrong.');
    end
    if ~isVector && size(third_input) ~= m 
        error('Normalization vector has a wrong size.') 
    end
    
    % Assing the 3rd input
    if isVector
        zeroStrategy = third_input;
    else
        norm_factor = third_input;
    end
end



%% Main

if ~isVector
    % Relative Error (MATRIX)
    err_rel = zeros(n, 1);
    for i = 1 : n
        err_rel(i,1) = norm( Y(:,i) - Y_full(:,i) ) / norm_factor(i);
    end

else
    % Relative Error (VECTOR)
    err_rel = abs(Y - Y_full) ./ (eps + Y_full);
    
    % Zero Strategy
    I = find(Y_full < 3 * eps);
    if zeroStrategy
        err_rel(I) = abs(Y(I) - Y_full(I));
    end
    
end



end







