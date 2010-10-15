% 
% Usage:   alpha=mexNesterov(X,D,param);
%
% Name: mexNesterov
%
% WARNING: This function has not been tested intensively
%
% Description: mexNesterov is an implementation of 
%     two algorithms presented in 
%     "Gradient methods for minimizing composite objective function" by
%     Y. Nesterov, 2007
%     
%     It is designed for solving a large number of small or medium-sized
%     decomposition problem (and not for a single large one).  It first
%     computes the Gram matrix D'D and then processes the  input signals in
%     parallel.  It aims at addressing the following problem
%
%     for all columns x_i of X, 
%        min_{alpha_i} 0.5||x_i-Dalpha_i||_2^2 + lambda1 ||alpha_i||_1 + ...
%             ... + 0.5lambda2||alpha_i||_2^2 + lambda3 FL(alpha_i)
%
%      When lambda3 == 0, it solves the Elastic Net.
%      When lambda3 == 0 and lambda2 == 0, it solves the Lasso
%
% Inputs: X:  double m x n matrix   (input signals)
%               m is the signal size
%               n is the number of signals to decompose
%         D:  double m x p matrix   (dictionary)
%               p is the number of elements in the dictionary
%         param: struct
%            param.lambda1  (parameter)
%            param.lambda2  (optional, parameter, 0 by default)
%            param.lambda3  (optional, parameter, 0 by default)
%            param.accelerated (optional, false by default) implements 
%              accelerating scheme of Nesterov's method.
%            param.rho  : Lipshitz constant estimate of the cost function.
%              It has not to be accurate. A line search takes care of 
%              tuning internally automatically this parameter.
%            param.tol :  (optional, 0.001 by default) tolerance parameter
%              for the stopping criterion
%            param.itermax (optional, 100 by default) maximum number 
%              of iterations
%            param.numThreads (optional, number of threads for exploiting
%              multi-core / multi-cpus. By default, it takes the value -1,
%              which automatically selects all the available CPUs/cores).
%
% Output: alpha: double p x n matrix (output coefficients)
%
% Note: this function admits a few experimental usages, which have not
%     been extensively tested:
%         - single precision setting (even though the output alpha is double 
%           precision)
%
% Author: Julien Mairal, 2009


