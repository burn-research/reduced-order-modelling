function xp = buildMap(xq, p, varargin)
%% Description
%{
 Build map or indeces for the x-p (spatial axis, parameter) space. 

 INPUTS: xq - spatial coordinates (CELL or DOUBLE);
         p  - parameters (DOUBLE, CELL);
%}


%% Checks

if (~isa(p,'double') && ~isa(p,'cell')) || (~isa(xq,'double') && ~isa(xq,'cell'))
    error('Inputs must be either a DOUBLE or CELL variable.');
end


%% Main

% Conditions
conds = [isa(p,'double'), isa(xq,'double')];

% Choose
if conds(1) && conds(2)         % True, True
    xp = Fun1(xq, p);
elseif ~conds(1) && conds(2)    % False, True
    xp = Fun2(xq, p);
elseif conds(1) && ~conds(2)    % True, False
    xp = Fun3(xq, p);
else                            % False, False
    xp = Fun4(xq, p);
end




end


% Both are DOUBLE
function xp = Fun1(xq, p, varargin)

% Sizes
[rows_p, cols_p] = size(p);      
[rows_xq, cols_xq] = size(xq);
 
% Initialize
xp1 = repmat(xq, rows_p, 1);
xp2 = zeros(rows_xq * rows_p, cols_p);

% Build the output
for i = 1 : rows_p
    temp = repmat(p(i,:), rows_xq, 1);  % Repeat this set of parameters
    
    % Get the starting and final rows to update
    i_start = 1 + (i-1) * rows_xq;
    i_final = i * rows_xq;
    
    xp2(i_start:i_final, :) = temp;     % Update
end

% Output
xp = [xp1, xp2]; 

end


% P is CELL, XQ is DOUBLE 
function xp = Fun2(xq, p, varargin)

% Sizes
[rows_p, cols_p] = size(p);      
[rows_xq, cols_xq] = size(xq);

repxq = repmat(xq, rows_p, 1);

xp1 = mat2cell(repxq, ones(rows_p * rows_xq,1)); 
xp2 = {};

for i = 1 : rows_p
    new = cellstr(p{i}); new = repmat(new, rows_xq, 1);
    xp2 = [xp2; new];
end

xp = [xp1, xp2];

end


% P is DOUBLE, XQ is CELL
function xp = Fun3(xq, p, varargin)

% Sizes
[rows_p, cols_p] = size(p);      
l = length(xq);
% Initialize
xp = [];
for j = 1 : length(xq)
    for i = 1 : l
        [rows, cols] = size( xq{i} );    
        temp = [xq{i}, repmat( p(i,:), rows, 1 )];
        xp = [xp; temp];
    end
end

end



% Both are CELL
function xp = Fun4(xq, p, varargin)



end

