function [w, beta, sigma2, logl] = reg_online(Y, Z, w_0, beta_0, sigma2_0, gam);

% reg_online     The on-line EM algorithm.
% Use: [w, beta, sigma2, logl] = reg_batch(Y, Z, w_0, beta_0, sigma2_0, gam) where
%       Y: responses (1 x n)
%       Z: covariates (d x n)
%       gam: step size (1 x n)
% and
%       w: (1 x m x n + 1)
%       beta: (d x m x n + 1)
%       sigma2 (1 x m x n + 1)
%
% $Id: reg_online.m,v 1.3 2008/03/17 11:13:18 cappe Exp $

% Parameters
INHIBIT_LAG = 20;
% In inhibition stage we let gamma(i) = 1/i and do not perfom update
gam(1:INHIBIT_LAG) = 1./(1:INHIBIT_LAG);

% Dimension
n = length(Y);
m = length(w_0);
d = size(Z,1);

% Initialization
w = zeros(1,m,n+1);
beta = zeros(d,m,n+1);
sigma2 = zeros(1,m,n+1);
w(:,:,1) = w_0;
beta(:,:,1) = beta_0;
sigma2(:,:,1) = sigma2_0;

logl = zeros(1,n);
bS1 = zeros(1,m);
bS2 = zeros(d,m);
bS3 = zeros(d,d,m);
bS4 = zeros(1,m);
for i = 1:n
    % The E-SA step
    [bS1_one, bS2_one, bS3_one, bS4_one, logl(i)] = ...
        reg_e(Y(i), Z(:,i), w(:,:,i), beta(:,:,i), sigma2(:,:,i));
    bS1 = (1-gam(i))*bS1 + gam(i)*bS1_one;
    bS2 = (1-gam(i))*bS2 + gam(i)*bS2_one;
    bS3 = (1-gam(i))*bS3 + gam(i)*bS3_one;
    bS4 = (1-gam(i))*bS4 + gam(i)*bS4_one;
    % M-step
    if (i >= INHIBIT_LAG)
      [w(:,:,i+1), beta(:,:,i+1), sigma2(:,:,i+1)] = reg_m(bS1, bS2, bS3, bS4);
    else
      % Inhibition phase
      w(:,:,i+1) = w(:,:,i);
      beta(:,:,i+1) = beta(:,:,i);
      sigma2(:,:,i+1) = sigma2(:,:,i);
    end
end
