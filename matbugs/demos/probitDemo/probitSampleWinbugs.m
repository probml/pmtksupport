function [samples, stats] = probitSampleWinbugs(X, y, Sigma, Nsamples, mu, useLogistic)
% Draw samples from beta ~ p(beta|X,y) 
% where p(beta) = N(beta,mu,Sigma)
% and p(y|X(i,:)=x, beta) = Bernoulli(Phi(beta'*x))
% where Phi is the cumulative Gaussian.
% If logistic=1, we use logit instead.
% We use Winbugs.
% Probit does not work well - often the chains don't mix.


[n p] = size(X);
if nargin < 5,  mu = zeros(p,1); end
if nargin < 6, useLogistic = 0;  end


dataStruct = struct('p', p, 'n', n, 'X', X, 'y', y, ...
                   'mu', mu, 'precMat', inv(Sigma));

Nchains = 5;

% we initialize the params randomly - probably a bad idea
clear initStructs
for i=1:Nchains
  %S.beta = mvnrnd(mu, Sigma);
  S.beta = mvnrnd(mu, eye(p));
  initStructs(i) = S;
end

if useLogistic
  fname = which('logitRegressionModel.txt');
else
  fname = which('probitRegressionModel.txt');
end

[samples, stats] = matbugs(dataStruct, fname, ...
		'init', initStructs, ...
		'view', 1, 'nburnin', 10000, 'nsamples', Nsamples, ...
		'thin', 1, 'nChains', Nchains, ...
		'monitorParams', {'beta'}, ...
		'Bugdir', 'C:/Program Files/WinBUGS14');

    
end