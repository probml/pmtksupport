% demo of finite mixture model
dd = size(wpcounts,1);
NN = size(wpcounts,2);
KK = 50;
aa = 2;
bb = 10000;
numiter = 100;

hh.dd = dd;
hh.aa = bb;

% construct data
xx = num2cell(wpcounts,1);

% initialize component assignment
zz = ceil(rand(1,NN)*KK);

% initialize finite mixture
dpm = dpm_init(KK,aa,Multinomial(hh),xx,zz);

% initialize records
record.KK = zeros(1,numiter);
record.zz = zeros(NN,numiter);

% run
tic; lasttime = toc;
for iter = 1:numiter
  if toc>lasttime+1, 
    lasttime = toc;
    fprintf(1,'DP mixture: iter# %d\n',iter);
    dpm_demonips_summarize(dpm)
  end

  % gibbs iteration 
  dpm = dpm_gibbs(dpm,1);

  % record keeping
  record.KK(iter) = sum(dpm.nn>0);
  record.zz(:,iter) = dpm.zz;
end

