clear all

% Analyse finney's data - see homework 6 of
%http://www.stat.rutgers.edu/~madigan/bayes06/

V = [3.7,3.5,1.25,0.75,0.8,0.7,0.6,1.1,0.9,0.9,0.8,0.55,0.6,1.4,0.75,2.3,3.2,...
0.85,1.7,1.8,0.4,0.95,1.35,1.5,1.6,0.6,1.8,0.95,1.9,1.6,2.7,2.35,1.1,1.1,1.2,0.8,...
0.95,0.75,1.3];
R = [0.825,1.09,2.5,1.5,3.2,3.5,0.75,1.7,0.75,0.45,0.57,2.75,3.0,...
2.33,3.75,1.64,1.6,1.415,1.06,1.8,2.0,1.36,1.35,1.36,1.78,1.5,1.5,1.9,0.95,0.4,...
0.75,0.03,1.83,2.2,2.0,3.33,1.9,1.9,1.625];
y = [1,1,1,1,1,1,0,0,0,0,0,0,0,1,1,1,1,1,0,1,0,0,0,0,1,0,1,0,1,0,1,0,0,1,1,1,0,0,1];

x1 = log(10*V);
x2 = log(10*R);

n = length(x1)
X = [ones(n,1) x1(:) x2(:)];

p = 3;
Sigma = 1000*eye(p);
mu = zeros(p,1);
     
useLogistic = 1;
[samples, stats] = probitSampleWinbugs(X, y, Sigma, Nsamples, mu, useLogistic);

Nchains = size(samples.beta, 1);

figure(1); clf
colors = 'rgbkcy';
for c=1:Nchains
  for i=1:min(p, 9)
    subplot(3,3,i)
    plot(samples.beta(c,:,i), colors(c));
    hold on
  end
end
title(sprintf('beta'))

%% average across chains
%%samplesBugs = squeeze(mean(samples.beta, 1))'; % samples(1:p, 1:S)

% concatenate chains
keepChains = 3; % in case some haven't mixed
keepChains = 1:Nchains;
samplesBugs = reshape(samples.beta(keepChains,:,:), ...
		      length(keepChains)*Nsamples, p)'; %samples(1:p, 1:N*S)

nr = 2; nc = 2;

figure(2); clf
for i=1:min(p, nr*nc)
  subplot(nr, nc,i)
  ksdensity(samplesBugs(i, :))
  title(sprintf('mean = %5.3f, sd = %5.3f', ...
		mean(samplesBugs(i,:)), std(samplesBugs(i,:))))
end
%set(gcf, 'name', 'winbugs')




if 0
figure(3);
for i=1:Nsamples
  mse(i) = norm(beta - mean(samplesBugs(:, 1:i),2));
end
plot(mse)
title('winbugs (beta - E[beta|Nsamples])^2')
end
