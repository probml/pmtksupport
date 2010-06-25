function [bS1, bS2, bS3, bS4, logl] = reg_e(Y, Z, w, beta, sigma2);

% reg_e         Computes the EM statistics for one data point.
% Use: [bS1, bS2, bS3, bS4, logl] = reg_e(Y, Z, w, beta, sigma2) where
%       Y:      response (scalar)
%       Z:      covariates (d x 1)
%       w:      weights (1 x m)
%       beta:   regression coeffs. (d x m)
%       sigma2: variances (1 x m)
% and
%       bS1 (1 x m)
%       bS2 (d x m)
%       bS3 (d x d x m)
%       bS4 (1 x m)
%
% $Id: reg_e.m,v 1.1 2007/09/04 16:37:10 cappe Exp $

% Constants
SQRT_2PI = 2.506628274631;

% Dimensions
d = length(Z);
m = length(w);

% Unormalized weights
bw = exp( -0.5 * ( Y(ones(1,m)) - Z'*beta ).^2 ./ sigma2 ) .* w ./ sqrt(sigma2) / SQRT_2PI;
% Log-likelihood
lkhd = sum(bw);
logl = log(lkhd);
% Normalized weights
bw = bw / lkhd;

% Statistics
bS1 = bw;
bS2 = Y*Z*bw;
bS3 = zeros(d,d,m);
bS3(:,:,1) = Z*Z';
for j = 2:m
  bS3(:,:,j) = bw(j)*bS3(:,:,1);
end
bS3(:,:,1) = bw(1)*bS3(:,:,1);
bS4 = Y^2*bw;
