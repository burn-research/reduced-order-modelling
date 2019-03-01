function [gamma, debug] = cpca(X0, my_constraint, scores, modes, x_ave, x_sig, cpca_guess_corrector, options)
%% Description
% X0: (n_obs x n_var) data from which MODES and SCORES were extracted, which PCA was
%                    performed on.
% my_constraint: function handle of the non linear constraint to be
%                respected.
% scores: (n_obs x n_pc) PCA scores
% modes: PCA modes
% x_ave: vector of mean values
% x_sig: vector of scaling factors
%

%% Input
if isempty(my_constraint)
    my_constraint = @cpca_con;
end
if nargin < 7 || isempty(cpca_guess_corrector)
    cpca_guess_corrector = 1;
end
if nargin < 8
    options = [];
end

%% Evaluating the CPCA scores:
% Starting points (or guesses)
gamma0 = cpca_guess_corrector * scores; 
disp('--- Evaluating CPCA scores ---');
% Solver:
if isa(my_constraint, 'function_handle') 
    % Non-linear constrained optimization problem
    gamma = solv(my_constraint, scores', modes, [], [], gamma0', X0', x_ave, x_sig, options);
elseif isa(my_constraint, 'cell') || isa(my_constraint, 'double')
    % Linear constrained problem
    [gamma, debug] = solv2(obj, gamma0);
end

%% Output:
if nargout > 1
    varargout{1} = debug;
end
end


%------- Solving for the gamma's ---------------%
function [gamma] = solv(constraint_fun, scores, modes, lb, ub, gamma0, X, x_ave, x_sig, cpca_options)
n_obs = size(X, 2);
k = size(scores, 1);
% Evaluation of each gamma(:,i)
gamma = zeros(k, n_obs); % Memory allocation
% Possible algorithms for the solver of the constrained non-linear problem:
x = {'interior-point', 'sqp'};  % {'interior-point'; 'sqp'; 'active-set'};
% Loop repetitions
parfor i = 1 : n_obs       
    % Starting guess:
    guess = gamma0(:,i);
    % Solving
    j = 1;
    goforth = true;
    while goforth && (j <= length(x))
        % Options:
        options = setOptions(x{j});
        % Use "FMINCON":
        [gammasol, fval, exitflag] = fmincon(...
            @(gamma) ObjFun(gamma, X(:,i), modes),... % Objective function
            guess,...       % Starting guess
            [],[],...       % Linear inequality costraints
            [],[],...       % Linear equality costraints 
            lb,ub,...       % Lower and upper bounds (lb, ub)
            @(gamma) constraint_fun(gamma, modes, x_ave, x_sig, cpca_options),... % Non-linear function for costraints
            options);       % Options
        guess = gammasol;
        % Exitflag
        if exitflag < 1
            j = j + 1;
        else
            goforth = false;
        end 
    end
    % Solution:
    gamma(:,i) = gammasol; 
    % Show real time to prove work is in progress
    fprintf('Solved for case %i out of %i\n', i, n_obs);
end % parfor
end % end of solv()
%------- End: Solving for the gamma's ---------------% 


%------- Set options ---------------%
function options = setOptions(s)
switch s
    case 'interior-point'
        MaxIter = 5e5;     TolFun = 1e-12;     TolX = 1e-12; 
        TolCon = 1e-8; 
    case 'sqp'
        MaxIter = 1e6;     TolFun = 1e-15;     TolX = 1e-15;
        TolCon = 1e-12;
    case 'active-set'
        MaxIter = 5e2;     TolFun = 1e-18;     TolX = 1e-18;
        TolCon = 1e-12;
end
MaxFunEvals = 1e11;
options = optimoptions('fmincon',...
        'Algorithm',s,...
        'MaxFunEvals',MaxFunEvals, 'MaxIter',MaxIter, 'TolFun',TolFun,...
        'TolX',TolX, 'TolCon',TolCon,...
        'Display','off'); % 'iter', 'final-detailed'
end
%------- End: Set options ---------------%


%------- ObjFun ----------------------%
function f = ObjFun(gamma, y, modes) 
y_rec = modes * gamma;
err = y - y_rec;
h = norm(y) + eps;
f = .5 * norm( err ) / ( h );
end
%------- End: ObjFun ---------------%

%------- cpca_con ------------------%
function [c, ceq] = cpca_con(cpca_scores, modes, x_ave, x_sig, options)
%% Description:
% This is an non-linear constraints which requires the rebuilt data vector
% y to have all positive values. It can be a common constraint to be used,
% so it is present here by default.
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
% cpca_scores: vector of CPCA scores
% modes: PCA modes
% x_ave: mean values (vector)
% x_sig: scaling factors (vector)
%

%% Inputs
if nargin < 3 || isempty(x_ave)
    x_ave = 0;
end
if nargin < 4 || isempty(x_sig)
    x_sig = 1;
end

%% Pre-proc:
a_tol = 1e-16;
x_ave(x_ave < a_tol) = 0; % To avoid problems
% Column vectors
cpca_scores = cpca_scores(:);
x_ave = x_ave(:);
x_sig = x_sig(:);

%% Non-linear inequalities:
% The vector y is rebuilt
y = modes * cpca_scores;
if isempty(options)
    y = x_ave + x_sig .* y;
else
    if isfield(options, 'idx') && isfield(options, 'k') && ~isempty(options.idx) && ~isempty(options.k)
        y = scale_center_decode2(y', options.mesh, options.variables, 1, options.sc, options.scale_on_Q, options.idx, options.k)';
    else
        y = scale_center_decode(y', options.mesh, options.variables, 1, options.sc, options.scale_on_Q)';
    end
end
% All values have to be > 0
c = - min(y);                       

%% Non-linear equalities:
ceq = []; % No non-linear inequality constraints.

end
%------- cpca_con ------------------%


function [Y0, sc] = scale_center_main(Y, mesh, variables, X, cent_crit, scal_crit, scale_on_Q)
sc.Y_ave = [];
sc.Y_gamma = [];
sc.Q_ave = [];
sc.Q_gamma = [];
if scale_on_Q == 1
    Q = getStateVars(Y', mesh, variables, X)';
    [Q0, sc.Q_ave] = center(Q, cent_crit); % Y_ave is a vector
    [Q0, sc.Q_gamma] = scale(Q0, Q, scal_crit); % Y_gamma is a vector
    Y0 = leaveStateVars(Q0', mesh, variables, X)'; clear Q; clear Q0
elseif scale_on_Q == 0
    [Y0, sc.Y_ave] = center(Y, cent_crit); % Y_ave is a vector
    [Y0, sc.Y_gamma] = scale(Y0, Y, scal_crit); % Y_gamma is a vector
elseif scale_on_Q == 2
    [Y0, sc.Y_ave] = center(Y, cent_crit);
    Q = getStateVars(Y', mesh, variables, X)';
    Q0 = getStateVars(Y0', mesh, variables, X)';
    [Q0, sc.Q_gamma] = scale(Q0, Q, scal_crit);
    Y0 = leaveStateVars(Q0', mesh, variables, X)'; clear Q; clear Q0
end
end
function [Y_out] = scale_center_decode(Y_rec, mesh, variables, X, sc, scale_on_Q)
if scale_on_Q == 1
    Q = getStateVars(Y_rec', mesh, variables, X)';
    Q = unscale(Q, sc.Q_gamma);
    Q = uncenter(Q, sc.Q_ave);
    Y_out = leaveStateVars(Q', mesh, variables, X)';
elseif scale_on_Q == 0
    Y_out = unscale(Y_rec, sc.Y_gamma);
    Y_out = uncenter(Y_out, sc.Y_ave);
elseif scale_on_Q == 2
    Q = getStateVars(Y_rec', mesh, variables, X)';
    Q = unscale(Q, sc.Q_gamma);
    Y_out = leaveStateVars(Q', mesh, variables, X)';
    Y_out = uncenter(Y_out, sc.Y_ave);
end
end
function [Y_out] = scale_center_decode2(Y_rec, mesh, variables, X, sc, scale_on_Q, idx, k)
I = (idx == k);
if scale_on_Q == 1
    Q = getStateVars(Y_rec', mesh, variables, X)';
    Q = unscale(Q, sc.Q_gamma);
    Q = uncenter(Q, sc.Q_ave);
    Y_out = leaveStateVars(Q', mesh, variables, X)';
elseif scale_on_Q == 0
    Y_out = unscale(Y_rec, sc.Y_gamma(I));
    Y_out = uncenter(Y_out, sc.Y_ave(I));
elseif scale_on_Q == 2
    Q = getStateVars(Y_rec', mesh, variables, X)';
    Q = unscale(Q, sc.Q_gamma);
    Y_out = leaveStateVars(Q', mesh, variables, X)';
    Y_out = uncenter(Y_out, sc.Y_ave(I));
end
end



