function [mu, lambda, nu] = bayes_linear_post(X, w, V, a_tau, b_tau)
%% [mu, lambda, nu] = bayes_linear_post(X, w, V, a_tau, b_tau)
%
% returns the posterior for bayes_linear_fit(_ard), given the inputs x being
% the rows of X.
%
% The arguments are the ones returned by bayes_linear_fit(_ard), specifying
% the parameter posterior
%
% N(w1 | w, tau^-1 V) Gam(tau | a_tau, b_tau).
%
% The predictive posteriors are of the form
%
% St(y | mu, lambda, nu),
%
% which is a Student's t distribution with mean mu, precision lambda, and nu
% degrees of freedom. All of mu and lambda a vectors, one per input x. nu
% is a scalar as it is the same for all x.

mu = X * w;
lambda = (a_tau / b_tau) ./ (1 + sum(X .* (X * V), 2));
nu = 2 * a_tau;
