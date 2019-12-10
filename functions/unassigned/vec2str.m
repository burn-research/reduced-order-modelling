function str = vec2str(vec, varargin)

n_args = length(varargin);

[rows, cols] = size(vec);
if rows ~= 1 && cols ~= 1
    error('Works only for vectors');
end

if rows == 1 && cols == 1
    str = num2str(vec);
    return
end

sep = '_';
if n_args > 0
    sep = varargin{1};
end

str = [];
ll = length(vec);
for ii = 1 : ll
    str = [num2str(vec(ii)), sep, str];
end

end


