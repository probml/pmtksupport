% 
% Usage:   alpha=mexOMP(X,D,param);
%
% Name: mexOMP
%
% Description: mexOMP is an efficient implementation of the
%     Orthogonal Matching Pursuit algorithm. It is optimized
%     for solving a large number of small or medium-sized 
%     decomposition problem (and not for a single large one).
%     It first computes the Gram matrix D'D and then perform
%     a Cholesky-based OMP of the input signals in parallel.
%     It aims at addressing the following NP-hard problem
%     
%     for all columns x_i of X, 
%         min_{alpha_i} ||alpha_i||_0  s.t  ||x_i-Dalpha_i||_2^2 <= eps
%         or
%         min_{alpha_i} ||x_i-Dalpha_i||_2^2  s.t. ||alpha_i||_0 <= L
%         
%
% Inputs: X:  double m x n matrix   (input signals)
%            m is the signal size
%            n is the number of signals to decompose
%         D:  double m x p matrix   (dictionary)
%            p is the number of elements in the dictionary
%            All the columns of D should have unit-norm !
%         param: struct
%            param.L  (maximum number of elements in each decomposition)
%            param.eps (threshold on the squared l2-norm of the residual
%            param.numThreads (optional, number of threads for exploiting
%            multi-core / multi-cpus. By default, it takes the value -1,
%            which automatically selects all the available CPUs/cores).
%
% Output: alpha: double sparse p x n matrix (output coefficients)
%
% Note: this function admits a few experimental usages, which have not
%     been extensively tested:
%      - single precision setting (even though the output alpha is double 
%        precision)
%      - Passing an int32 vector of length n to param.L allows to provide
%        a different parameter L for each input signal x_i
%      - Passing a double vector of length n to param.eps allows to provide
%        a different parameter eps for each input signal x_i
%
% Author: Julien Mairal, 2009


