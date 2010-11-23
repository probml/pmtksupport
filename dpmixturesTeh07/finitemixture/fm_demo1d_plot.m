function pp = fm_demo1d_plot(fm,yy);
% plots data and density, assumes in 1D
% evaluates density at yy

tolerance = .001;
ih = ishold;

xx = cat(2,fm.xx{:});
plot(xx,zeros(size(xx)),'x');
hold on

% draw mixing proportions over active and inactive components
mpi = dirrnd([fm.nn+fm.aa/fm.KK]);

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

plot(yy,pp);

if ~ih, hold off; end
