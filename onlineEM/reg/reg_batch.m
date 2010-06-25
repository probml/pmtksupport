function [w, beta, sigma2, logl] = reg_batch(Y, Z, w_0, beta_0, sigma2_0, nit);

% reg_batch     The usual (batch) EM algorithm.
% Use: [w, beta, sigma2, logl] = reg_batch(Y, Z, w_0, beta_0, sigma2_0, nit) where
%       Y: responses (1 x n)
%       Z: covariates (d x n)
% and
%       w: (1 x m x nit+1)
%       beta: (d x m x nit+1)
%       sigma2 (1 x m x nit+1)
% Note that the last iteration is here just to compute the likelihood.
%
% $Id: reg_batch.m,v 1.1 2007/09/04 16:37:10 cappe Exp $

% Dimension
n = length(Y);
m = length(w_0);
d = size(Z,1);

% Initialization
w = zeros(1,m,nit+1);
beta = zeros(d,m,nit+1);
sigma2 = zeros(1,m,nit+1);
w(:,:,1) = w_0;
beta(:,:,1) = beta_0;
sigma2(:,:,1) = sigma2_0;

logl = zeros(1,nit+1);
for iit = 1:nit+1
  bS1 = zeros(1,m);
  bS2 = zeros(d,m);
  bS3 = zeros(d,d,m);
  bS4 = zeros(1,m);
  % Accumulate statistics over the E-step
  for i = 1:n
    [bS1_one, bS2_one, bS3_one, bS4_one, logl_one] = ...
        reg_e(Y(i), Z(:,i), w(:,:,iit), beta(:,:,iit), sigma2(:,:,iit));
    bS1 = bS1 + bS1_one;
    bS2 = bS2 + bS2_one;
    bS3 = bS3 + bS3_one;
    bS4 = bS4 + bS4_one;
    logl(iit) = logl(iit) + logl_one;
  end
  % Normalize
  bS1 = bS1/n;
  bS2 = bS2/n;
  bS3 = bS3/n;
  bS4 = bS4/n;
  % M-step
  if (iit <= nit)
    [w(:,:,iit+1), beta(:,:,iit+1), sigma2(:,:,iit+1)] = reg_m(bS1, bS2, bS3, bS4);
  end
end

% Also normalize likelihoods
logl = logl/n;
