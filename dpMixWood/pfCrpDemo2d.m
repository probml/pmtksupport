% Use particle filtering to fit a DP mixture of Gaussians 
% to 2d Old Faithful data

setSeed(0);
X = loadData('faithful'); % 272x2
X = standardizeCols(X);
%X = X(1:150,:);

% The model requires us to specify:
%
%   \Lambda_0^{-1} : our prior for per-neuron spike variability about it's
%                    mean spike
%   v_0 : how much confidence we have in this prior (must be greater than
%   spike_dimensions_to_retain for mathematical reasons)
%   mu_0 : the mean waveform for a particular cell
%   k_0 : how much confidence we have in this prior
%   a_0, b_0 : gamma(a_0, b_0) on the CRP concentration parameter alpha
%

D = 2;
mu_0 = zeros(D,1);
k_0 = .01;
lambda_0 = 0.75*eye(D);
v_0 = (D+1);
a_0 = 1;
b_0 = 1;

if 0
% sample from prior 
for trial=1:3
[samples labels means covariances] = sample_igmm_prior(50,a_0,b_0,mu_0,lambda_0,k_0,v_0);
scatter(samples(:,1),samples(:,2),[],labels);
drawnow
end
end

alpha = 2;
num_particles = 200;
partial_labels = 1;
[z, w, K] = particle_filter(X',num_particles, a_0, b_0, mu_0, k_0, v_0, ...
    lambda_0, alpha, partial_labels);

[Np Nd] = size(z);
Kstar = double(max(K));
post = zeros(Nd, Kstar);
for i=1:Nd
  post(i,:) = normalize(hist(double(z(:,i)), 1:Kstar));
end
[junk, yhat] = max(post,[],2);

if 0
figure; hold on
colors = pmtkColors;
for k=1:Kstar
  ndx = find(yhat==k);
  h=plot(X(ndx,1), X(ndx,2),'o');
  set(h, 'color', colors{k});
end
end

[styles, colors, symbols, str] =  plotColors;
Ns = [5 10 20];
for N=Ns(:)'
figure; hold on
for i=1:N
  [pmax k] = max(post(i,:));
  h=plot(X(i,1), X(i,2), symbols(k));
  set(h, 'color', colors(k));
  ms  = 10*pmax;
  ms = min(ms, 20);
  ms = max(ms, 5);
  set(h, 'markersize', ms); 
end
title(sprintf('N=%d', N))
printPmtkFigure(sprintf('pfCrpDemo2dN%d', N));
end

