function [wght, rate, wght_stats, rate_stats, logl] = ...
   pm_rem_step(count, wght, rate,  wght_stats, rate_stats, gamma, dnorm, update);

%pm_rem_step
%         Elementary step of the REM algorithm.
%         Use: [wght, rate, wght_stats, rate_stats, logl] = pm_rem_step ...
%           (count, wght, rate,  wght_stats, rate_stats, gamma, dnorm, update);
%         note: does not perform any checks.
%
% $Id: pm_rem_step.m,v 1.2 2006/04/04 17:05:39 oKp Exp $

%%% 1: Compute a posteriori probabilities and likelihood
% Compute all densities
postprob = exp(-rate + count*log(rate) - dnorm);
% Compute unormalized a posteriori probability
postprob = postprob .* wght;
% Compute log-likelihood (not that this is computed prior to updating the
% parameters)
logl = log(sum(postprob));
% Normalization
postprob = postprob / sum(postprob);

%%% 2: Update the statistics
wght_stats = gamma*postprob + (1-gamma)*wght_stats;
rate_stats = gamma*(postprob*count) + (1-gamma)*rate_stats;

if (update)
  %%% 3: Update the parameters (M step formula)
  wght = wght_stats;
  rate = rate_stats./wght_stats;
end
