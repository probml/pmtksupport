% 
% Usage:  U=mexSparseProject(B,param);
%
% Name: mexSparseProject
%
% Description: mexSparseProject solves various optimization 
%     problems, including projections on a few convex sets.
%     It aims at addressing the following problems
%     for all columns b of B in parallel
%       1) when param.mode=1 (projection on the l1-ball)
%           min_u ||b-u||_2^2  s.t.  ||u||_1 <= thrs
%       2) when param.mode=2
%           min_u ||b-u||_2^2  s.t. ||u||_2^2 + lambda1||u||_1 <= thrs
%       3) when param.mode=3
%           min_u ||b-u||_2^2  s.t  ||u||_1 + 0.5lambda1||u||_2^2 <= thrs 
%       4) when param.mode=4
%           min_u 0.5||b-u||_2^2 + lambda1||u||_1  s.t  ||u||_2^2 <= thrs
%       5) when param.mode=5
%           min_u 0.5||b-u||_2^2 + lambda1||u||_1 +lambda2 FL(u) + ... 
%                                                   0.5lambda_3 ||u||_2^2
%          where FL denotes a "fused lasso" regularization term.
%       6) when param.mode=6
%          min_u ||b-u||_2^2 s.t lambda1||u||_1 +lambda2 FL(u) + ...
%                                             0.5lambda3||u||_2^2 <= thrs
%       7) when param.mode=7
%           min_u ||b-u||_2^2  s.t  lambda_1||u||_1 + (1-lambda1)||u||_2^2 <= thrs 
%           
%        When param.pos=true and param.mode <= 4,
%        it solves the previous problems with positivity constraints 
%
% Inputs: B:  double m x n matrix   (input signals)
%               m is the signal size
%               n is the number of signals to project
%         param: struct
%           param.thrs (parameter)
%           param.lambda1 (parameter)
%           param.lambda2 (parameter)
%           param.lambda3 (parameter)
%           param.mode (see above)
%           param.pos (optional, false by default)
%           param.numThreads (optional, number of threads for exploiting
%             multi-core / multi-cpus. By default, it takes the value -1,
%             which automatically selects all the available CPUs/cores).
%
% Output: U: double m x n matrix (output matrix)
%
% Note: this function admits a few experimental usages, which have not
%     been extensively tested:
%         - single precision setting 
%
% Author: Julien Mairal, 2009


