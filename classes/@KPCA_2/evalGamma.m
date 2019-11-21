function [gamma, varargout] = evalGamma(obj, varargin)
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


%% Evaluating the gamma's:

% Starting points (or guesses)
guesses = obj.guess_corrector * obj.a; 


disp(' '); disp('%%% - Solving: evaluating gamma`s... - %%%');


% Solver:
temp = obj.mycon;
if isa(temp,'function_handle')
    % Non-linear constrained optimization problem
    [gamma, debug] = solv(obj, [], [], debug, guesses, fval_precisionTarget);
    
elseif isa(temp,'cell') || isa(temp,'double')
    % Linear constrained problem
    [gamma, debug] = solv2(obj, guesses);
    
end


%% Output:
    if (nargout>1)
        varargout{1} = debug;
    end

    
end %end of evalGamma






%------- Solving for the gamma's ---------------%
function [gamma, debug] = solv(obj, lb, ub, debug, guesses, varargin)


fval_precisionTarget = Inf;
if ~isempty(varargin)
    fval_precisionTarget = varargin{1};
end


% Load the needed PCA modes
modes = obj.modes(:, 1:obj.k); 


% Evaluation of each gamma(:,i)
gamma = zeros(obj.k, obj.cols); %Memory allocation


for i = 1 : obj.cols   
    
    % Possible algorithms for the solver of the constrained non-linear problem:
    x = {'interior-point', 'sqp'};  % {'interior-point'; 'sqp'; 'active-set'};
    
    % Target function:
    y = obj.Ycs(:,i);
    
    
    % Objective function's coefficients:
    objfun.y = y;
    objfun.modes = modes;
    objfun.sqmodes = modes'*modes;

    
    % Starting guess:
    guess = guesses(:,i);
    
    
    %% Solver:
    j = 1;
    
    
    % Constraint violation check [REMEMBER c(x) <= 0]:
    [c, ceq] = obj.mycon(obj.a(:,i), obj);
    
    
    % If the constraint is violated (or the guess is not the set of PCA 
    % scores), proceed with CPCA 
    if ( ~isempty(c) && any(c > 0) ) || obj.guess_corrector ~= 1
        goforth = true;
    else
        goforth = false;
        
        gammasol = obj.a(:,i); % Return this solution
        
        % Set the debug info for what happened
        outString = 'No fmincon was run';
        exitflag = outString; 
        fval = outString;
        output = outString; 
        lambda = outString; 
        grad = outString;
        hessian = outString; 
        
        fprintf(['\nNo fmincon was run for the sample ',num2str(i),'\n']);
    end
    
    
    % Solving
    tic;
    while goforth && j <= length(x) % Use x{3} only if things are easy
        
        % Options:
        options = setOptions(x{j});
        
        %always at least tries 'interior-point'
        disp(' '); 
        disp(['- Using: "',x{j},'" for the ',num2str(i),'-th training point -']);
        
        % Use "fmincon":
        [gammasol,fval,exitflag,output,lambda,grad,hessian] = fmincon(...
            @(gamma) ObjFun(gamma, objfun),... % Objective function
            guess,...       % Starting guess
            [],[],...       % Linear inequality costraints
            [],[],...       % Linear equality costraints 
            lb,ub,...       % Lower and upper bounds (lb, ub)
            @(gamma) obj.mycon(gamma,obj),... % Non-linear function for costraints
            options);       % Options
        guess = gammasol;
        
        disp(' '); disp(['Done using [',num2str(j),'] "',x{j},'".']);
        
        
        if exitflag < 1
            j = j + 1;
        else
            goforth = false;
        end
        
    end % while 
    
    
    % Time evaluation
    elapsed_time = toc;
    fprintf('\nElapsed time was %d s.\n', elapsed_time);
    
    
    % Solution:
    gamma(:,i) = gammasol;    
    
    
    % If FVAL is too large (if precision to recover the real data is low)
    if ~isa(grad,'char') && (norm(fval) > fval_precisionTarget)
        gamma(:,i) = obj.a(:,i);
        warning('CPCA scores had a high FVAL. Set equal to the PCA scores.');
    end
    
    
    %% Debug Variables
    debug.Exitflag{i} = exitflag; 
    debug.Fval{i} = fval;
    debug.Output{i} = output; 
    debug.Lambda{i} = lambda; 
    debug.Grad{i} = grad; 
    debug.Hessian{i} = hessian; 
    
    
    %% Display progress
    disp(['Progress: ', num2str(100*i/obj.cols), ' %']);
    
    
end % for

end % end of solv()

%------- End: Solving for the gamma's ---------------% 




%------- Set options ---------------%
function options = setOptions(s)

TolCon = 1e-12;

switch s
    case 'interior-point'
        MaxIter = 5e2;     TolFun = 1e-11;     TolX = 1e-11; 
        TolCon = 1e-7; 
    case 'sqp'
        MaxIter = 1e4;     TolFun = 1e-15;     TolX = 1e-15;
    case 'active-set'
        MaxIter = 5e2;     TolFun = 1e-18;     TolX = 1e-18; 
end

MaxFunEvals = 1e10;

options = optimoptions('fmincon',...
        'Algorithm',s,...
        'MaxFunEvals',MaxFunEvals, 'MaxIter',MaxIter, 'TolFun',TolFun,...
        'TolX',TolX, 'TolCon',TolCon,...
        'Display','final'); %iter, final-detailed

end
%------- End: Set options ---------------%




%------- ObjFun ----------------------%
function f = ObjFun(gamma, objfun) 
%% Remember



%% Evaluation
    y = objfun.y;
    y_rec = objfun.modes * gamma;
    
    err = y - y_rec;
    h = norm(objfun.y) + eps;
    
    f = .5 * norm( err ) / ( h );
    
end
%------- End: ObjFun ---------------%



%     objfun.c3 = y'*y; objfun.c2 = -2*(y'*modes); objfun.c1 = 1;




