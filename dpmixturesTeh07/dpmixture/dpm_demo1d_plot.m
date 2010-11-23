function pp = dpm_plot1d(dpm,yy);
% plots data and density, assumes in 1D
% evaluates density at yy

tolerance = .001;
ih = ishold;

xx = cat(2,dpm.xx{:});
plot(xx,zeros(size(xx)),'x');
hold on

% draw mixing proportions over active and inactive components
mpi = dirrnd([dpm.nn dpm.aa]);

% draw mixing proportions over inactive components by stick-breaking.
while (mpi(end)>tolerance)
  bb = betarnd(1,dpm.aa);
  mpi([end end+1]) = mpi(end)*[bb 1-bb];
end

% draw mean and variance for each component
KK = length(mpi);
mu = zeros(1,KK);
sigma = zeros(1,KK);
for kk = 1:KK
  ll = min(kk,dpm.KK);
  [mu(kk) sigma(kk)] = rand(dpm.qq{min(kk,dpm.KK+1)});
end

% compute density
pp = zeros(size(yy)); 
for kk = 1:KK
  pp = pp + mpi(kk)*1/sqrt(2*pi*sigma(kk))*exp(-.5*(yy-mu(kk)).^2/sigma(kk));
end

plot(yy,pp);

if ~ih, hold off; end
