% generate DP by repeatedly splitting pi in 2.
% cheats by drawing from full Dirichlet first

ll = 10;
nn = 2^15;
aa = 3;

gg = gamrnd(aa/nn,1,1,nn);
dd = gg/sum(gg);
xx = (0:nn)/nn*ll;
yy = dd;

for ii = 1:log2(nn)
  zz{ii} = dd;
  dd = sum(reshape(dd,[2 length(dd)/2]),1);
end
zz{log2(nn)+1} = dd;
zz = zz(end:-1:1);

aa = [0-.1 10+.1 0-.1 40];
clf
title('A Draw from a Dirichlet process');
for ii = 2.^(0:log2(nn))
  xx = (0:ii)/ii*ll;
  yy = zz{log2(ii)+1}*ii/ll;
  cc = [.7 1 1]*(1-(log2(ii)+1)/(log2(nn)+1))/1.5;
  hh = patch([0 xx(floor((1:2*ii)/2)+1) ll],[0 yy(ceil((1:2*ii)/2)) 0],cc);
  pause
  set(hh,'edgecolor',cc);
end

