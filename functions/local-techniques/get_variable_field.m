function [field, idx] = get_variable_field(Y, mesh, list_of_var, var, X, x, varargin)
%% Description:
% Y: obs x var, snapshot matrix
%

%% Input
if ~exist('X', 'var') || isempty(X)
    X = [];
    x = [];
end

%% Main
if ischar(var)
    [field, idx] = get_variable_field_local(Y, mesh, list_of_var, var, X, x);
elseif iscell(var)
    field = [];
    idx = [];
    len = length(var);
    for ii = 1 : len
        [temp, tmp] = get_variable_field_local(Y, mesh, list_of_var, var{ii}, X, x);
        field = [field, temp];
        idx = [idx; tmp];
    end
elseif isnumeric(var)
    field = [];
    idx = [];
    len = length(var);
    for ii = 1 : len
        [temp, tmp] = get_variable_field_local(Y, mesh, list_of_var, list_of_var{var(ii)}, X, x);
        field = [field, temp];
        idx = [idx; tmp];
    end
else
    error('There seems to be a problem with the input VARIABLE.');
end
end

function [field, idx] = get_variable_field_local(Y, mesh, vars, variable, X, x, varargin)
% Y: obs x var 
n_mesh = size(mesh,1);
n_vars = length(vars);
rows = size(Y, 1);
for ii = 1 : n_vars
    if strcmp(variable, vars{ii})
        jj = ii;
        break
    end
end
x1 = (jj - 1) * n_mesh + 1;
x2 = jj * n_mesh;
idx = x1 : x2;
field = Y(:, idx);
if isempty(X)
    return
end
I = ismember(X, x, 'rows');
field = field(I,:);
end



