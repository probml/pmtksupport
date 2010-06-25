function [R,S] = cov2cor(V)
% Compute the correlation matrix from the covariance matrix
% Sigma(i) = sqrt( Covariance(i,i) );
% Corr(i,j) = Covariance(i,j)/( Sigma(i)*Sigma(j) );

S = sqrt(diag(V));
R = V ./ (S*S');

