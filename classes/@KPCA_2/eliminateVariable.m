function eliminateVariable(obj, varargin)

% Number of inputs
n_args = length(varargin);

% Input checks
if n_args == 0
    error('Specify which variable to eliminate.');
end
var = varargin{1};

% Recursive if there is more than one variable to eliminate
if iscell(var)
    % In the recursive run, VAR is of type CHAR, this IF won't activate
    l = length(var);
    for i = 1 : l
        obj.eliminateVariable(var{i}); 
    end
    return
end

% Get variable index
var_I = get_var_pos(var, obj.variable_names);

% Find the rows to eliminate
n = size(obj.training_data, 1);
m = size(obj.mesh,1);
I = true(n,1);
if obj.is_mesh_variable
    I(1+m*(var_I-1):m*var_I) = false;
else
    I(var_I) = false;
end

% Update training data
obj.training_data = obj.training_data(I,:);
obj.variable_names(var_I) = [];
if ~isempty(obj.original_data)
    obj.original_data = obj.original_data(I,:);
end
end


function i = get_var_pos(var, vars)

i = 0;
for ii = 1 : length(vars)
    if strcmp(var, vars{ii})
        i = ii;
        return
    end
end

end