function [training_points, prediction_points, training_data, original_data]...
    = swapSample(points, mesh, imv, variable_names, ...
    training_points, prediction_points, training_data, varargin)
%% Description

%% Input 
n_points = size(points, 1); % Number of points
p = points(1,:);

n_args = length(varargin);
if n_args > 0
    original_data = varargin{1};
else
    original_data = [];
end
    
%% Find the point
mesh_dim = size(mesh, 2);
if imv
    x = training_points;
    x_test = prediction_points;
else
    x = unique( training_points(:,mesh_dim+1:end) );
    x_test = unique( prediction_points(:,mesh_dim+1:end) );
end
[flag_train, flag_pred] = ismember_here(x, x_test, p);

% If both are true or false, there is a problem
if any(flag_train) + any(flag_pred) ~= 1
    fprintf('\nPoint not found or found in both sets.\n');
    return
end

%% Swap
if any(flag_train)
    % Swap to prediction points
    if imv
        [training_points, prediction_points, training_data, original_data] = swap(...
            training_points, prediction_points, training_data, original_data, flag_train);
    else
        [data, data_test] = getYform(x, x_test, training_data, original_data, mesh, variable_names);
        [x, x_test, data, data_test] = swap(x, x_test, data, data_test, flag_train);
        [training_points, prediction_points, training_data, original_data] ...
            = getQform(x, x_test, data, data_test, mesh, variable_names);
    end
elseif any(flag_pred)
    % Swap to training points
    if imv
        [prediction_points, training_points, original_data, training_data] ...
            = swap(prediction_points, training_points, original_data, training_data, flag_pred);
    else
        [data, data_test] = getYform(x, x_test, training_data, original_data, mesh, variable_names);
        [x_test, x, data_test, data] = swap(x_test, x, data_test, data, flag_pred);
        [training_points, prediction_points, training_data, original_data] ...
            = getQform(x, x_test, data, data_test, mesh, variable_names);
    end
end

% Recursive call
if n_points > 1
    [training_points, prediction_points, training_data, original_data] = ...
        swapSample(points(2:end,:), mesh, imv, variable_names,... 
        training_points, prediction_points, training_data, original_data);
end

end

function [x, x_to, data, data_to] = swap(x, x_to, data, data_to, flag)

data_to = [data_to, data(:,flag)];
data(:,flag) = [];
x_to = [x_to; x(flag,:)];
x(flag,:) = [];
[~, I] = sortrows(x_to);
x_to = x_to(I,:);
data_to = data_to(:,I);

end 

function [Y, Y_test] = getYform(xp, xpk, Q, Q_test, mesh, vars)

Y = leaveStateVars(Q, mesh, vars, xp);
Y_test = leaveStateVars(Q_test, mesh, vars, xpk);

end

function [xp, xpk, Q, Q_test] = getQform(x, x_test, Y, Y_test, mesh, vars)

[Q, xp] = getStateVars(Y, mesh, vars, x);
[Q_test, xpk] = getStateVars(Y_test, mesh, vars, x_test);

end

function [flag_train, flag_pred] = ismember_here(x, x_test, p)

% If there are no test points, set them to Inf: this way they can't be
% selected as they will be too far away
if isempty(x_test)
    x_test = Inf * ones(size(x));
end

% Get distances
d_train = x - repmat(p, size(x,1), 1);
d_train = sqrt( sum(d_train.^2, 2) );
d_test = x_test - repmat(p, size(x_test,1), 1);
d_test = sqrt( sum(d_test.^2, 2) );

% Find where the smallest distance is
d_train = d_train(:);
d_test = d_test(:);
[~, pos] = min([d_train; d_test]);

flag_train = false(length(d_train),1);
flag_pred = false(length(d_test),1);
if pos <= length(d_train)
    flag_train(pos) = true;
else
    flag_pred(pos-length(d_train)) = true;
end

end


