clear all;
randn('seed',0);
% Data are generated
X=randn(64,1000);
D=randn(64,256);
D=D./repmat(sqrt(sum(D.^2)),[size(D,1) 1]);

% parameter of the optimization procedure are chosen
param.L=10; % not more than 10 non-zeros coefficients
param.eps=0.1; % squared norm of the residual should be less than 0.1
param.numThreads=-1; % number of processors/cores to use; the default choice is -1
                    % and uses all the cores of the machine

tic
alpha=mexOMP(X,D,param);
t=toc

fprintf('%f signals processed per second\n',size(X,2)/t);
