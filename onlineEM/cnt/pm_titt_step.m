function [wght_out, rate_out, logl] = ...
   pm_titt_step(count, wght_in, rate_in, gamma, dnorm, update);

%pm_titt_step
%         Elementary step of the recursive EM algorithm a la Titterington.
%         Use: [wght, rate, logl] = pm_rem_step ...
%           (count, wght, rate, gamma, dnorm, update);
%         note: does not perform any checks.
%
% $Id: pm_titt_step.m,v 1.1 2006/04/06 08:36:51 oKp Exp $

%%% 1: Compute a posteriori probabilities and likelihood
% Compute all densities
postprob = exp(-rate_in + count*log(rate_in) - dnorm);
% Compute unormalized a posteriori probability
postprob = postprob .* wght_in;
% Compute log-likelihood (not that this is computed prior to updating the
% parameters)
logl = log(sum(postprob));
% Normalization
postprob = postprob / sum(postprob);

if (update)
  %%% 2: Update the parameters
  wght_out = gamma*postprob + (1-gamma)*wght_in;
  % For the rates that care of possible negative values
  g_tmp = gamma*postprob./wght_in;
  rate_out = g_tmp*count + (1-g_tmp).*rate_in;
  if (any(rate_out < 0))
    fprintf('Warning: rejected update\n');
    rate_out = rate_in;
  end
end
