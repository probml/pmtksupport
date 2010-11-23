function pp = fm_demo1d_density(fm,yy);
% pp = fm_demo1d_density(fm,yy)
% evaluates density at yy on the current sample from the finite mixture.

% draw mixing proportions over active and inactive components
mpi = dirrnd(fm.nn+fm.aa/fm.KK);

% draw mean and variance for each component
KK = length(mpi);
mu = zeros(1,KK);
sigma = zeros(1,KK);
for kk = 1:KK
  [mu(kk) sigma(kk)] = rand(fm.qq{kk});
end

% compute density
pp = zeros(size(yy)); 
for kk = 1:KK
  pp = pp + mpi(kk)*1/sqrt(2*pi*sigma(kk))*exp(-.5*(yy-mu(kk)).^2/sigma(kk));
end
