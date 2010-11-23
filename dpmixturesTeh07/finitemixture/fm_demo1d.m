% demo of finite mixture model in 1D
dd = 1;
KK = 2;
NN = 10;
xx = [-2+.5*randn(1,5) 2+.5*randn(1,5)];
aa = 1;
s0 = 3;
ss = 1;
numiter = 1000;

hh.dd = dd;
hh.ss = s0^2/ss^2;
hh.vv = 5;
hh.VV = ss^2*eye(dd);
hh.uu = zeros(dd,1);


% set range over which to evaluate density
yy = -15:.01:15;

% initialize records
record.KK = zeros(1,numiter);
record.p0 = zeros(numiter,length(yy)); % densities
record.pp = zeros(numiter,length(yy)); % densities

% construct finite mixture with no data items to get prior samples of densities
fm0 = fm_init(2,aa,Gaussian(hh),{},[]);

% initialize finite mixture with data items.
xx = num2cell(xx); % data
zz = ceil(rand(1,NN)*KK); % initialize component assignment
fm = fm_init(KK,aa,Gaussian(hh),xx,zz);

% run
for iter = 1:numiter
  fprintf(1,'finite mixture: iter# %d\r',iter);
  drawnow; 

  % gibbs iteration 
  fm = fm_gibbs(fm,1);

  % record keeping
  record.KK(iter) = sum(fm.nn>0);
  record.pp(iter,:) = fm_demo1d_density(fm,yy); 
  record.p0(iter,:) = fm_demo1d_density(fm0,yy);
end
fprintf(1,'\n');

figure(1);
fm_demo1d_summarize(fm0,yy,record.p0);
title('Prior over densities');
figure(2);
fm_demo1d_summarize(fm,yy,record.pp);
title('Posterior over densities');
