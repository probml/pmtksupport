function [w, beta, sigma2] = reg_m(bS1, bS2, bS3, bS4);

% reg_m         The M-step
% Use: [w, beta, sigma2] = reg_m(bS1, bS2, bS3, bS4).
%
% $Id: reg_m.m,v 1.1 2007/09/04 16:37:10 cappe Exp $

% Dimensions
[d, m] = size(bS2);

w = bS1;
beta = zeros(d,m);
sigma2 = zeros(1,m);
for j = 1:m
  beta(:,j) = bS3(:,:,j) \ bS2(:,j);
  sigma2(j) = (bS4(j) - beta(:,j)'*bS2(:,j))/bS1(j);
end
