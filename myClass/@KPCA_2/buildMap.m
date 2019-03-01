function xp = buildMap(xq, p)
%% Description
%{
 Build map or indeces for the x-p (spatial axis, parameter) space
 INPUTS: xq - spatial coordinates;
         p  - parameter or another set of quantities;
%}


%%

if ~isa(p,'double') && ~isa(p,'cell')
    error('Second input must be either a DOUBLE or CELL variable.');
end


m = length(p);      % Number of parameter points
n = length(xq);     % Number of spatial points


if isa(p,'double')
    xp = doubleFun(xq, p, m, n);  % if p is a DOUBLE vector
elseif isa(p,'cell') 
    xp = cellFun(xq, p, m, n);    % if p is CELL array (of variable names)
end


end



function xp = doubleFun(xq, p, m, n)

% Initializing variables
xp1 = repmat(xq, m, 1);     xp2 = xp1 * 0; 

for i = 1 : m
    xp2(1+n*(i-1):(n*i), 1) = p(i) * ones(n,1);
end

% Indeces for the x-p space
xp = [xp1, xp2]; 

end



function xp = cellFun(xq, p, m, n)

repxq = repmat(xq, m, 1);

xp1 = mat2cell(repxq, ones(m*n,1)); 
xp2 = {};

for i = 1 : m
    new = cellstr(p{i}); new = repmat(new, n, 1);
    xp2 = [xp2; new];
end

xp = [xp1, xp2];

end


