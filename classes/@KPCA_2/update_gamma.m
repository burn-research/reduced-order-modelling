function update_gamma(obj, varargin)

% Check dataset was provided
n_args = length(varargin);
if n_args < 1
    error('The CPCA optimization problem requires a dataset that the CPCA formulation needs to approximate.\n');
end

% Evaluate the CPCA scores
fprintf('Evaluation of the CPCA scores is in progress.');
obj.cpca_scores = evalGamma(obj, varargin{1});
fprintf('Evaluation of the CPCA scores has terminated.');

end


function [gamma, varargout] = evalGamma(obj, min_data, varargin)
%% Input
fval_precisionTarget = Inf;
if ~isempty(varargin)
    fval_precisionTarget = varargin{1};
end

%% Debug structure:
debug.Fval = {}; 
debug.Output = {}; 
debug.Lambda = {}; 
debug.Grad = {}; 
debug.normGrad = {}; 
debug.Hessian = {}; 
debug.HesCheck = {};  

%% Evaluating the CPCA scores:
k = obj.pca_approximation_order;
% Starting points (or guesses)
guesses = obj.cpca_guess_corrector * obj.pca_scores(1:k,:); 

disp('--- Evaluating CPCA scores ---');
% Solver:
temp = obj.my_constraint;
if isa(temp,'function_handle') 
    % Non-linear constrained optimization problem
    [gamma, debug] = solv(obj, [], [], debug, guesses, min_data, fval_precisionTarget);
elseif isa(temp,'cell') || isa(temp,'double')
    % Linear constrained problem
    [gamma, debug] = solv2(obj, guesses);
end

%% Output:
if nargout > 1
    varargout{1} = debug;
end
    
end % end of evalGamma()


%------- Solving for the gamma's ---------------%
function [gamma, debug] = solv(obj, lb, ub, debug, guesses, min_data, varargin)

fval_precisionTarget = Inf;
if ~isempty(varargin)
    fval_precisionTarget = varargin{1};
end

% Load the needed PCA modes
modes = obj.pca_modes(:, 1:obj.pca_approximation_order); 

% Evaluation of each gamma(:,i)
gamma = zeros(obj.pca_approximation_order, size(obj.training_data, 2)); % Memory allocation

% Possible algorithms for the solver of the constrained non-linear problem:
x = {'interior-point', 'sqp'};  % {'interior-point'; 'sqp'; 'active-set'};

% Pointer to function
constraint_fun = @obj.my_constraint;

% PCA scores
scores = obj.pca_scores(1:obj.pca_approximation_order, :);

% Loop repetitions
n_obs = size(obj.training_data, 2);

parfor i = 1 : n_obs       
    % Starting guess:
    guess = guesses(:,i);
    
    % Constraint violation check [REMEMBER c(x) <= 0]:
    [c, ceq] = constraint_fun(scores(:,i), obj.mean_column, obj.scaling_factors, modes);
    
    % If the constraint is violated (or the guess is not the set of PCA 
    % scores), proceed with CPCA 
    if ( ~isempty(c) && any(c > 0) ) || obj.cpca_guess_corrector ~= 1
        goforth = true; % Continue
    else
        goforth = false; % Stop
        gammasol = scores(:,i); % Return this solution
        
        % Set the debug info for what happened
        outString = 'No fmincon was run';
        exitflag = outString; 
        fval = outString;
    end
    
    % Solving
    j = 1;
    while goforth && (j <= length(x))
        % Options:
        options = setOptions(x{j});
        
        % Use "FMINCON":
        [gammasol, fval, exitflag] = fmincon(...
            @(gamma) ObjFun(gamma, min_data(:,i), modes),... % Objective function
            guess,...       % Starting guess
            [],[],...       % Linear inequality costraints
            [],[],...       % Linear equality costraints 
            lb,ub,...       % Lower and upper bounds (lb, ub)
            @(gamma) constraint_fun(gamma, obj.mean_column, obj.scaling_factors, modes),... % Non-linear function for costraints
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

% Evaluation
y_rec = modes * gamma;

err = y - y_rec;
h = norm(y) + eps;

f = .5 * norm( err ) / ( h );
    
end
%------- End: ObjFun ---------------%


