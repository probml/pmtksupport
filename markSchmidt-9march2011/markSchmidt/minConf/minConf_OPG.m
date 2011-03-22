function [x,f,funEvals,projects] = minConf_OPG(funObj,x,funProj,options)
% function [x,f] = minConF_OPG(funObj,x,funProj,options)
%
% Function for using Optimal Projected Gradient to solve problems of the form
%   min funObj(x) s.t. x in C
%
%   @funObj(x): function to minimize (returns gradient as second argument)
%   @funProj(x): function that returns projection of x onto C
%
%   options:
%       verbose: level of verbosity (0: no output, 1: final, 2: iter (default), 3:
%       debug)
%       optTol: tolerance used to check for progress (default: 1e-6)
%       maxIter: maximum number of calls to funObj (default: 500)
%       numDiff: compute derivatives numerically (0: use user-supplied
%       derivatives (default), 1: use finite differences, 2: use complex
%       differentials)
%       suffDec: sufficient decrease parameter in Armijo condition (default
%       : 1e-4)
%       interp: type of interpolation (0: step-size halving, 1: quadratic,
%       2: cubic)
%       memory: number of steps to look back in non-monotone Armijo
%       condition
%       curvilinear: backtrack along projection Arc (default: 0)
%       testOpt: test optimality condition (default: 1)
%       feasibleInit: if 1, then the initial point is assumed to be
%       feasible
%
%   Notes: 
%       - if the projection is expensive to compute, you can reduce the
%           number of projections by setting testOpt to 0


nVars = length(x);

% Set Parameters
if nargin < 4
    options = [];
end
[verbose,numDiff,optTol,progTol,maxIter,suffDec,interp,feasibleInit,testOpt,L] = ...
    myProcessOptions(...
    options,'verbose',2,'numDiff',0,'optTol',1e-5,'progTol',1e-9,'maxIter',500,'suffDec',1e-4,...
    'interp',2,'feasibleInit',0,...
    'testOpt',1,'L',[]);

% Output Log
if verbose >= 2
    if testOpt
        fprintf('%10s %10s %10s %15s %15s %15s\n','Iteration','FunEvals','Projections','Step Length','Function Val','Opt Cond');
    else
        fprintf('%10s %10s %10s %15s %15s\n','Iteration','FunEvals','Projections','Step Length','Function Val');
    end
end

% Make objective function (if using numerical derivatives)
funEvalMultiplier = 1;
if numDiff
    if numDiff == 2
        useComplex = 1;
    else
        useComplex = 0;
    end
    funObj = @(x)autoGrad(x,useComplex,funObj);
    funEvalMultiplier = nVars+1-useComplex;
end

% Evaluate Initial Point
if ~feasibleInit
    x = funProj(x);
end
[f,g] = funObj(x);
projects = 1;
funEvals = 1;

% Optionally check optimality
if testOpt
    projects = projects+1;
    if max(abs(funProj(x-g)-x)) < optTol
        if verbose >= 1
        fprintf('First-Order Optimality Conditions Below optTol at Initial Point\n');
        end
        return;
    end
end

% Initialize
mu = 0;
gamma = 1;
alphap = 1;

if isempty(L)
    t = min(1,1/sum(abs(g)));
    L = 1/t;
end
y = x;
f_y = f;
g_y = g;

i = 1;
while funEvals <= maxIter
    
    while 1
        b = -gamma+mu;
        alpha= (b+ sqrt(b*b + 4* L * gamma)) / (2*L);
        beta= (gamma - gamma* alphap) / (alphap * gamma + alphap* L * alpha);
        
        if i > 1
            y = x + beta*(x - x_old);
            [f_y,g_y] = funObj(y);
            funEvals = funEvals+1;
        end
        
        x_new = funProj(y - g_y/L);
        projects = projects+1;
        [f_new,g_new] = funObj(x_new);
        funEvals = funEvals+1;
        
        % Backtrack if not below Lipschitz
        
        l_sum = f_new - f_y - g_y'*(x_new-y);
        r_sum = (1/2)*(x_new-y)'*(x_new-y);
        if l_sum <= r_sum*L
            break;
        else
            if verbose == 3
                fprintf('Decreasing Step Size\n');
            end
            
            if isLegal(l_sum)
                L = max(2*L,l_sum/r_sum);
            else
                L = 2*L;
            end
            
            if max(abs(g_y/L)) < progTol
                if verbose == 3
                    fprintf('Line search failed\n');
                end
                return
            end
        end
    end
    
    % Take step
    f_old = f;
    g_old = g;
    x_old = x;
    x = x_new;
    f = f_new;
    g = g_new;
    
    if testOpt
        optCond = max(abs(funProj(x-g)-x));
        projects = projects+1;
    end

    % Output Log
    if verbose >= 2
        if testOpt
            fprintf('%10d %10d %10d %15.5e %15.5e %15.5e\n',i,funEvals*funEvalMultiplier,projects,1/L,f,optCond);
        else
            fprintf('%10d %10d %10d %15.5e %15.5e\n',i,funEvals*funEvalMultiplier,projects,1/L,f);
        end
    end

    % Check optimality
    if testOpt
        if optCond < optTol
            if verbose >= 1
            fprintf('First-Order Optimality Conditions Below optTol\n');
            end
            break;
        end
    end

    if sum(abs(x-x_old)) < progTol
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

    if funEvals*funEvalMultiplier > maxIter
        if verbose >= 1
            fprintf('Function Evaluations exceeds maxIter\n');
        end
        break;
    end
    
    
    gamma=L* alpha* alpha;    alphap=alpha;
    ratio=L / (l_sum/ r_sum);
    
    if (ratio > 5)
        if verbose >= 3
        fprintf('Increasing step size\n');
        end
        L=L*0.8;
    end
    if ratio < 1
        fprintf('WTF!\n');
        %pause
    end
    
    i = i + 1;
end