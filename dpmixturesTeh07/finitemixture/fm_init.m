function fm = fm_init(KK,aa,q0,xx,zz);
% initialize finite mixture model, with 
% KK mixture components,
% aa concentration parameter,
% q0 an empty component with hh prior,
% xx data, x_i=xx{i}
% zz initial cluster assignments (between 1 and KK).

fm.KK = KK;
fm.NN = length(xx);
fm.aa = aa;
fm.qq = cell(1,KK);
fm.xx = xx;
fm.zz = zz;
fm.nn = zeros(1,KK);

% initialize mixture components
for kk = 1:KK,
  fm.qq{kk} = q0;
end

% add data items into mixture components
for ii = 1:fm.NN
  kk = zz(ii);
  fm.qq{kk} = additem(fm.qq{kk},xx{ii});
  fm.nn(kk) = fm.nn(kk) + 1;
end

