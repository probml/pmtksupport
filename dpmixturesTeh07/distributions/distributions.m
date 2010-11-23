%
% Interface for Exponential Family Mixture Components with Conjugate Priors
% -------------------------------------------------------------------------
%
% This section of the provided code implements mixture components belonging to
% an exponential family with conjugate priors on the parameters.  The components
% are implemented as classes in matlab, with the following methods:
% 
% q = Gaussian(h) or q = Multinomial(h)
%         Class constructor for the two classes.  help Gaussian or
%         Multinomial for more information.
% q = additem(q,x)
% q = delitem(q,x)
%         Adds or deletes data item x from the component q.
% l = logmarginal(q)
%         Log marginal probability of data items in component q.
% l = logpredictive(q,x)
%         Log predictive probability of data item x given other data items in
%         component q.
% n = numitems(q)
%         Number of data items in component q.
% [mu,sigma] = Gaussian/map(q)
%         Returns the MAP estimate of the mean and covariance of data items.
% pi = Multinomial/mean(q)
%         Returns the mean estimate of the multinomial probabilities given
%         data items.
