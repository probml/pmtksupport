function dpm = dpm_init(KK,aa,q0,xx,zz);
% initialize DP mixture model, with 
% KK active mixture components,
% aa concentration parameter,
% q0 empty component with hh prior,
% xx data, x_i=xx{i}
% zz initial cluster assignments (between 1 and KK).

dpm.KK = KK;
dpm.NN = length(xx);
dpm.aa = aa;
dpm.qq = cell(1,KK+1);
dpm.xx = xx;
dpm.zz = zz;
dpm.nn = zeros(1,KK);

% initialize mixture components
% component KK+1 takes care of all inactive components
for kk = 1:KK+1,
  dpm.qq{kk} = q0;
end

% add data items into mixture components
for ii = 1:dpm.NN
  kk = zz(ii);
  dpm.qq{kk} = additem(dpm.qq{kk},xx{ii});
  dpm.nn(kk) = dpm.nn(kk) + 1;
end

