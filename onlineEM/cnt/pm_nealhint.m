function [wght, rate, logl] = pm_nealhint(count, wght, rate, Nit, alpha);

%pm_nealhint Estimates the parameters of a Poisson mixture using the incremenatl EM algorithm.
%         Use: [wght,rate,logl] = pm_nealhint(count,wght_0,rate_0,Nit,alpha).

% H2M Toolbox, Version 1.5
% Olivier Cappé, 30/12/97 - 31/12/97
% ENST Dpt. TSI / LTCI (CNRS URA 820), Paris

% Check input arguments
error(nargchk(4, 5, nargin));
if (nargin < 5)
  alpha = 1;
end
% Data length
T = length(count);
if (any(count < 0) | any(count ~= fix(count)))
  error('Data does not contain positive integers.');
end
count = reshape(count, T, 1);
% Number of mixture components
N = pm_chk(wght, rate);
wght = reshape(wght, 1, N);
rate = reshape(rate, 1, N);

% Variables
logl = zeros(1,Nit);
postprob = zeros(T, N);
% Compute log(count!), the second solution is usually much faster
% except if max(count) is very large
cm = max(count);
if (cm > 50000)
  dnorm = gammaln(count + 1);
else
  tmp = cumsum([0; log((1:max(count)).')]);
  dnorm = tmp(count+1);
end

% Main loop of the EM algorithm
for nit = 1:Nit
  if (nit == 1)
    % The first EM iteration is a usual
    %%% 1: Compute a posteriori probabilities and likelihood (E)
    % Compute all densities
    postprob = exp(-ones(T,1)*rate + count*log(rate) - dnorm*ones(1,N));
    % Compute unormalized a posteriori probability
    postprob = postprob .* (ones(T,1) * wght);
    % Compute loglikelihood
    logl(nit) = sum(log(sum(postprob')));
    fprintf(1, 'Iteration %d:\t%.3f\n', (nit-1), logl(nit));
    % Normalization
    postprob = postprob ./ ((sum(postprob'))' * ones(1,N));

    %%% 2: Reestimating mixture parameters (M)
    % Unormalized weights
    wght = sum(postprob);
    % Intensities
    rate = (count.' * postprob) ./ wght;
    % Normalize weigths
    wght = wght / sum(wght);
  else
    for nobs = 1:T
      % Replace just one element
      postprob(nobs,:) = exp(-rate + count(nobs)*log(rate) - dnorm(nobs));
      % Compute unormalized a posteriori probability
      postprob(nobs,:) = postprob(nobs,:) .* wght;
      % Normalization
      postprob(nobs,:) = postprob(nobs,:) / sum(postprob(nobs,:));
      if (alpha > 1)
        postprob(nobs,:) = alpha*postprob(nobs,:);
        postprob([1:(nobs-1), (nobs+1):T],:) = (1-(alpha-1)/(T-1)) ...
          *  postprob([1:(nobs-1), (nobs+1):T],:);
      end
      %%% 2: Reestimating mixture parameters (M) - THIS ONE IS VERY INEFICIENT
      % Unormalized weights
      wght = sum(postprob);
      % Intensities
      rate = (count.' * postprob) ./ wght;
      % Normalize weigths
      wght = wght / sum(wght);
    end
    % Compute all densities - JUST for LL
    postprob_logl = exp(-ones(T,1)*rate + count*log(rate) - dnorm*ones(1,N));
    % Compute unormalized a posteriori probability
    postprob_logl = postprob_logl .* (ones(T,1) * wght);
    % Compute loglikelihood
    logl(nit) = sum(log(sum(postprob_logl')));
    fprintf(1, 'Iteration %d:\t%.3f\n', (nit-1), logl(nit));
  end
end
