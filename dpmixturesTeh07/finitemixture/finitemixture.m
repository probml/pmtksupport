% Code to implement finite mixtures.  Assumes exponential family likelihood
% with conjugate prior over parameters.  Integrates out mixing proportions,
% and component parameters, Gibbs sample cluster indicator variables.
%
% Main files you can run:
% fm_demo1d   - 1D Gaussian data demo, plotting prior and posterior densities
% fm_demo2d   - 2D Gaussian data demo, plotting Gaussians as ellipses 
% fm_demonips - clustering NIPS papers demo
%
% Main file:
% fm_gibbs    - implements collapsed Gibbs sampling for finite mixture model.
% fm_init     - initializes a finite mixture model.
%
% Variables in finite mixture model:
% KK - number of active clusters
% NN - number of data items
% aa - alpha parameter of Dirichlet prior over mixing proportions
% qq - row cell vector of mixture components (MATLAB objects in distributions/)
% xx - row cell vector of data items
% zz - row vector of cluster indicator variables
% nn - row vector of number of data items per cluster
%
% Function interfaces:
% fm = fm_init(KK,aa,q0,xx,zz);
% fm = fm_gibbs(fm,numiter);

