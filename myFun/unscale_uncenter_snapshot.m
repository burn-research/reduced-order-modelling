function [Y_out] = unscale_uncenter_snapshot(Y, xq, variables, X, sc, scale_on_Q)
%% Description:
% INPUT
% Y: data matrix, obs x vars

%% Main
m = size(X, 1);
if scale_on_Q == 1
    Q = getStateVars(Y', xq, variables, X)';
    Q = unscale(Q, sc.Q_gamma);
    Q = uncenter(Q, sc.Q_ave);
    Y_out = getSnapshotMatrix(Q, xq, m)';
elseif scale_on_Q == 0
    Y_out = unscale(Y, sc.Y_gamma);
    Y_out = uncenter(Y_out, sc.Y_ave);
elseif scale_on_Q == 2
    Q = getStateVars(Y', xq, variables, X)';
    Q = unscale(Q, sc.Q_gamma);
    Y_out = getSnapshotMatrix(Q, xq, m)';
    Y_out = uncenter(Y_out, sc.Y_ave);
end
end