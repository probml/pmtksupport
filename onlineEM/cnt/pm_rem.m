function [wght, rate, logl] = pm_rem(count, wght_0, rate_0, gamma, update)

%pm_rem	  Estimates the parameters of a Poisson mixture using the REM algorithm.
%         Use: [wght,rate,logl] = pm_rem(count,wght_0,rate_0,gamma,update)
%         where wght and rate are the estimated model parameters, logl
%         contains the log-likehood values for the successive
%         iterations.
%
% $Id: pm_rem.m,v 1.3 2006/04/06 08:40:10 oKp Exp $

% Check input arguments
error(nargchk(4, 5, nargin));
% Data length
T = length(count);
if (any(count < 0) | any(count ~= fix(count)))
  error('Data does not contain positive integers.');
end
count = reshape(count, T, 1);
if (nargin < 5)
  update = ones(T,1);
  % The parameters should never be updated after the first step, otherwise the
  % rates will collapse to an identical value
  update(1) = 0;
end
% Number of mixture components
N = pm_chk(wght_0, rate_0);
wght_0 = reshape(wght_0, 1, N);
rate_0 = reshape(rate_0, 1, N);
% Stepsizes
if (length(gamma) ~= T)
  error('Incorrect number of stepsizes.');
end
if (any(gamma < 0) | any(gamma > 1))
  error('Stepsizes must be between 0 and 1');
end

% Output variables
wght = zeros(T,N);
rate = zeros(T,N);
logl = zeros(T,1);
% Compute log(count!), the second solution is usually much faster
% except if max(count) is very large
cm = max(count);
if (cm > 50000)
  dnorm = gammaln(count + 1);
else
  tmp = cumsum([0; log((1:max(count)).')]);
  dnorm = tmp(count+1);
end

% Main loop of the REM algorithm
for t = 1:T
  if (t == 1)
    wght_stats = zeros(1,N);
    rate_stats = zeros(1,N);    
    [wght(1,:), rate(1,:), wght_stats, rate_stats, logl(1)] = ...
      pm_rem_step(count(1), wght_0, rate_0, wght_stats, rate_stats, ...
      gamma(1), dnorm(1), update(1));
  else
    [wght(t,:), rate(t,:), wght_stats, rate_stats, logl(t)] = ...
      pm_rem_step(count(t), wght(t-1,:), rate(t-1,:), wght_stats, ...
      rate_stats, gamma(t), dnorm(t), update(t));   
  end
end
