function [c, ceq] = allPositiveCon(cpca_scores, mean_vec, scaling_factors, modes, varargin)
%% Description:
% This is an non-linear constraints which requires the rebuilt data vector
% y to have all positive values. It can be a common constraint to be used,
% so it is present in this class folder by default.
%
% OUTPUTS:
% c(x) is the array of nonlinear inequality constraints at x: 
% fmincon attempts to satisfy c(x) <= 0 for all entries of c.
%
% ceq(x) is the array of nonlinear equality constraints at x. 
% fmincon attempts to satisfy ceq(x) = 0 for all entries of ceq.
% 
% For example, x = fmincon(@myfun, x0, A, b, Aeq, beq, lb, ub, @mycon)
%
% INPUTS:
% gamma: vector of PCA scores
% m:     mean column vector
% d:     scaling factors
% modes: pca modes
%

%% Pre-proc:
a_tol = 1e-16;
mean_vec(mean_vec < a_tol) = 0; % To avoid problems

%% Additional inputs
nin = length(varargin);
if nin > 0 
    mean_vec = varargin{1};     % User-provided mean
    if nin > 1
        scaling_factors = varargin{2}; % User-provided scaling factors
    end
end

%% Non-linear inequalities:
% Relavant dimension
N = size(mean_vec, 1); 

% Diagonal matrix of scaling factors
D = spdiags(scaling_factors, 0, N, N);

% Unscale
modes = D * modes;

% The vector y is rebuilt
y = mean_vec + modes * cpca_scores; 

% All values have to be > 0
c = - min(y);                       

%% Non-linear equalities:
ceq = [];               % No non-linear inequality constraints.

end


