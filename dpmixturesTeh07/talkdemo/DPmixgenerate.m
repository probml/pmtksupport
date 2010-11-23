% generate from a dirichlet process mixture model

clear
aa = 1; % alpha
nn = 1000; % number of data points
ss = [1:10 50 100 200 500 1000]; % snapshot times
sigma = .5*eye(2); % mean covariance matrix
vsigma = 1;
dof = 10; % degree of freedom
mu  = zeros(2,1); % mean of means
mv  = 8*ones(2,1); % std of means
ax = 30;

T = [];
zz = zeros(2,nn);
for ii = 1:nn
  pp = [T aa];
  kk = sum(rand(1)*sum(pp) > cumsum(pp))+1;
  if kk < length(T)
    T(kk) = T(kk) + 1;
  else
    T(kk) = 1;
  end
  zz(ii) = kk;
end  

mm = zeros(2,length(T));
vv = zeros(2,2,length(T));
for kk = 1:length(T)
  mm(:,kk) = randn(2,1) .* mv + mu;
  vv(:,:,kk) = sqrtm(wishrnd(sigma,dof)) * sqrt(gamrnd(vsigma,1));
end

xx = zeros(2,nn);
for ii = 1:nn
  kk = zz(ii);
  xx(:,ii) = vv(:,:,kk) * randn(2,1) + mm(:,kk);
end

bb = 0:.02:2*pi;
figure(2)
for jj = 1:length(ss)
  sj = ss(jj);
  hh = hist(zz(1:sj),1:max(zz(1:sj)));
  cc = find(hh>=1);
  hh(jj) = plot(xx(1,1:sj),xx(2,1:sj),'.','markersize',7);
  hold on;
  for kk = cc
    uu = vv(:,:,kk);
    circ = mm(:,kk)*ones(1,length(bb)) + uu*[sin(bb);cos(bb)];
    plot(circ(1,:),circ(2,:),'linewidth',2,'color','k')
  end
  plot([-ax ax ax -ax -ax],[-ax -ax ax ax -ax]);
  hold off
  axis([-ax ax -ax ax])
  axis off
  pause
end
clear    


