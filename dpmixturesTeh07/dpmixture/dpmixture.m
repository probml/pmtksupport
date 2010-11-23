% Code to implement DP mixture models.  Assumes exponential family likelihood
% with conjugate prior over parameters.  Integrates out mixing proportions,
% and component parameters, Gibbs sample cluster indicator variables.
%
% Main files you can run:
% dpm_demo1d   - 1D Gaussian data demo, plotting prior and posterior densities
% dpm_demo2d   - 2D Gaussian data demo, plotting Gaussians as ellipses 
% dpm_demonips - clustering NIPS papers demo
%
% Main file:
% dpm_gibbs    - implements collapsed Gibbs sampling for DP mixture model.
% dpm_init     - initializes a DP mixture model.
%
% Variables in DP mixture model:
% KK - number of active clusters
% NN - number of data items
% aa - alpha parameter
% qq - row cell vector of mixture components (MATLAB objects in distributions/)
% xx - row cell vector of data items
% zz - row vector of cluster indicator variables
% nn - row vector of number of data items per cluster
%
% Function interfaces:
% dpm = dpm_init(KK,aa,q0,xx,zz);
% dpm = dpm_gibbs(dpm,numiter);

