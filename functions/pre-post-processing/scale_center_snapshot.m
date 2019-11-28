function [Y0, sc] = scale_center_snapshot(Y, xq, variables, X, cent_crit, scal_crit, scale_on_Q, mean_val, my_gamma)
%% Description
% Y: obs x vars

%% Inputs
if ~exist('mean_val', 'var') || isempty(mean_val)
    mean_val = [];
end
if ~exist('my_gamma', 'var') || isempty(my_gamma)
    my_gamma = [];
end
%% Main
sc.Y_ave = [];
sc.Y_gamma = [];
sc.Q_ave = [];
sc.Q_gamma = [];
m = size(X, 1);
if scale_on_Q == 1
    Q = getStateVars(Y', xq, variables, X)';
    [Q0, sc.Q_ave] = center(Q, cent_crit, mean_val); % Y_ave is a vector
    [Q0, sc.Q_gamma] = scale(Q0, Q, scal_crit, my_gamma); % Y_gamma is a vector
    Y0 = getSnapshotMatrix(Q0, xq, m)';
    clear Q; clear Q0
elseif scale_on_Q == 0
    [Y0, sc.Y_ave] = center(Y, cent_crit, mean_val); % Y_ave is a vector
    [Y0, sc.Y_gamma] = scale(Y0, Y, scal_crit, my_gamma); % Y_gamma is a vector
elseif scale_on_Q == 2
    [Y0, sc.Y_ave] = center(Y, cent_crit, mean_val);
    Q = getStateVars(Y', xq, variables, X)';
    Q0 = getStateVars(Y0', xq, variables, X)';
    [Q0, sc.Q_gamma] = scale(Q0, Q, scal_crit, my_gamma);
    Y0 = getSnapshotMatrix(Q0, xq, m)'; 
    clear Q; clear Q0
end
end
