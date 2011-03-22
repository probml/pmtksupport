function [x,f,funEvals] = minConf_QNST(funObj1,funObj2,x,funProj,options)

nVars = length(x);

if nargin < 5
    options = [];
end

[verbose,numDiff,optTol,progTol,maxIter,maxProject,suffDec,corrections,adjustStep,bbInit,...
    BBSToptTol,BBSTprogTol,BBSTiters,BBSTtestOpt] = ...
    myProcessOptions(...
    options,'verbose',2,'numDiff',0,'optTol',1e-5,'progTol',1e-9,'maxIter',500,'maxProject',100000,'suffDec',1e-4,...
    'corrections',10,'adjustStep',0,'bbInit',0,'BBSToptTol',1e-6,'BBSTprogTol',1e-10,'BBSTiters',10,'BBSTtestOpt',0);

% Output Parameter Settings
if verbose >= 3
   fprintf('Running QNST...\n');
   fprintf('Number of L-BFGS Corrections to store: %d\n',corrections);
   fprintf('Spectral initialization of BBST: %d\n',bbInit);
   fprintf('Maximum number of BBST iterations: %d\n',BBSTiters);
   fprintf('BBST optimality tolerance: %.2e\n',BBSToptTol);
   fprintf('BBST progress tolerance: %.2e\n',BBSTprogTol);
   fprintf('PQN optimality tolerance: %.2e\n',optTol);
   fprintf('PQN progress tolerance: %.2e\n',progTol);
   fprintf('Quadratic initialization of line search: %d\n',adjustStep);
   fprintf('Maximum number of function evaluations: %d\n',maxIter);
   fprintf('Maximum number of projections: %d\n',maxProject);
end

if verbose
    fprintf('%10s %10s %10s %15s %15s\n','Iteration','FunEvals','Projections','Step Length','Function Val');
end

% Evaluate Initial Objective
[f1,g] = funObj1(x);
f = f1+funObj2(x);
funEvals = 1;
projects = 0;

% Check optimality
optCond = max(abs(x-funProj(x-g,1)));
projects = 1;
if optCond < optTol
    if verbose >= 1
        fprintf('First-Order Optimality Conditions Below optTol at Initial Point\n');
    end
    return;
end

i = 1;
while 1
    
    if 0 % BBST
        if i == 1
            alpha = 1;
        else
            y = g-g_old;
            s = x-x_old;
            alpha = (s'*s)/(s'*y);
            if alpha <= 1e-10 || alpha > 1e10
                alpha = min(1,1/sum(abs(g)));
            end
        end
        p = funProj(x-alpha*g,alpha);
        projects = projects+1;
        
    else % QNST
        if i == 1
            p = funProj(x-g,1);
            projects = projects+1;
            S = zeros(nVars,0);
            Y = zeros(nVars,0);
            Hdiag = 1;
        else
            y = g-g_old;
        s = x-x_old;
        [S,Y,Hdiag] = lbfgsUpdate(y,s,corrections,verbose==3,S,Y,Hdiag);

        % Make Compact Representation
        k = size(Y,2);
        L = zeros(k);
        for j = 1:k
            L(j+1:k,j) = S(:,j+1:k)'*Y(:,j);
        end
        N = [S/Hdiag Y];
        M = [S'*S/Hdiag L;L' -diag(diag(S'*Y))];
        HvFunc = @(v)lbfgsHvFunc2(v,Hdiag,N,M);
        
        if bbInit
            % Use Barzilai-Borwein step to initialize sub-problem
            alpha = (s'*s)/(s'*y);
            if alpha <= 1e-10 || alpha > 1e10
                alpha = 1/norm(g);
            end
            
            % Solve Sub-problem
            xSubInit = funProj(x-alpha*g,alpha);
            projects = projects+1;
        else
            xSubInit = x;
        end
        % Solve Sub-problem
        [p,subProjects] = solveSubProblem(x,g,HvFunc,funObj2,funProj,BBSToptTol,BBSTprogTol,BBSTiters,BBSTtestOpt,xSubInit);
        projects = projects+subProjects;
        end
    end
    
   d = p-x;
   g_old = g;
   x_old = x;
    
   % Bound Step length on first iteration
   t = 1;
   if i == 1
      t = min(1,1/sum(abs(g)));
   end
   
   if t == 1
       x_new = p;
   else
    x_new = x+t*d;
   end
    [f1_new,g_new] = funObj1(x_new);
    f_new = f1_new + funObj2(x_new);
    funEvals = funEvals+1;

    f_old = f;
    while f_new > f
        if verbose
            fprintf('Backtracking\n');
        end
        t = t/2;
        x_new = x+t*d;
        [f1_new,g_new] = funObj1(x_new);
        f_new = f1_new + funObj2(x_new);
        funEvals = funEvals+1;
    end
    x = x_new;
    f = f_new;
    g = g_new;

        % Check Optimality
    optCond = max(abs(x-funProj(x-g,1)));
    projects = projects+1;
    
    if verbose
        fprintf('%10d %10d %10d %15.5e %15.5e %15.5e\n',i,funEvals,projects,t,f,optCond);
    end
    
    if optCond < optTol
        if verbose
            fprintf('First-order optimality below optTol\n');
        end
        break;
    end
    
    if max(abs(x-x_old)) < progTol
        if verbose >= 1
            fprintf('Step size below progTol\n');
        end
        break;
    end

    if abs(f-f_old) < progTol
        if verbose >= 1
            fprintf('Function value changing by less than progTol\n');
        end
        break;
    end

    if funEvals > maxIter
        if verbose
        fprintf('Exceeded maxIter funEvals\n');
        end
        break
    end

    i = i + 1;
end
end

function [p,subProjects] = solveSubProblem(x,g,H,funObj2,funProj,optTol,progTol,maxIter,testOpt,x_init)
% Uses BBST to solve for quasi-Newton soft-threshold direction
options.verbose = 0;
options.optTol = optTol;
options.progTol = progTol;
options.maxIter = maxIter;
options.testOpt = testOpt;

funObj = @(p)subHv(p,x,g,H);
[p,f,funEvals,subProjects] = minConf_BBST(funObj,funObj2,x_init,funProj,options);
end

function [f,g] = subHv(p,x,g,HvFunc)
d = p-x;
Hd = HvFunc(d);
f = g'*d + (1/2)*d'*Hd;
g = g + Hd;
end