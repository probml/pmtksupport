function pp = dpm_demo1d_density(dpm,yy);
% pp = dpm_post1d(dpm,yy)
% evaluates density at yy on the current sample from the DP mixture posterior.

% draw mixing proportions over active and inactive components
mpi = dirrnd([dpm.nn dpm.aa]);

% draw mixing proportions over inactive components by stick-breaking.
tolerance = .001;
while (mpi(end)>tolerance)
  bb = betarnd(1,dpm.aa);
  mpi([end end+1]) = mpi(end)*[bb 1-bb];
end

% draw mean and variance for each component
KK = length(mpi);
mu = zeros(1,KK);
sigma = zeros(1,KK);
for kk = 1:KK
  ll = min(kk,dpm.KK+1); % if kk>KK+1, is inactive component.
  [mu(kk) sigma(kk)] = rand(dpm.qq{ll});
end

% compute density
pp = zeros(size(yy)); 
for kk = 1:KK
  pp = pp + mpi(kk)*1/sqrt(2*pi*sigma(kk))*exp(-.5*(yy-mu(kk)).^2/sigma(kk));
end
