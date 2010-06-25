function [inf1, inf2] = reg_info(Y, Z, w, beta, sigma2);

% reg_info       Computes the observed information matrix for one data point (on
%                regression parameter only)
% Use: [inf1, inf2] = reg_info(Y, Z, w, beta, sigma2) where
%       Y:      response (scalar)
%       Z:      covariates (d x 1)
%       w:      weights (1 x m)
%       beta:   regression coeffs. (d x m)
%       sigma2: variances (1 x m)
%
% $Id: reg_info.m,v 1.1 2007/09/21 16:39:01 cappe Exp $

% Dimensions
m = length(w);
d = length(Z);

% E-step
post = reg_e(Y, Z, w, beta, sigma2);

M = Z*Z';
r = (Y(ones(1,m),:)-beta'*Z).^2;
% Actual observed information
inf1 = zeros(d,d,m);
for i = 1:m
  inf1(:,:,i) = M/sigma2(i)*post(i)*(1-r(i)/sigma2(i)*(1-post(i)));
end
% "Observation" at actual parameter value
inf2 = zeros(d,d,m);
for i = 1:m
  inf2(:,:,i) = M*r(i)/(sigma2(i)^2)*(post(i)^2);
end
