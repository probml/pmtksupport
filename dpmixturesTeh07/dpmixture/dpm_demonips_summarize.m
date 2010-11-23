function dpm_demonips_summasize(dpm);
% summarize the clusters in NIPS cluster demo.

% show clusters
for kk=1:dpm.KK
  fprintf(1,'cluster %d\n',kk);
  ii = find(dpm.zz==kk);
  for ll = 1:length(ii)
    fprintf(1,'%s\n',papers{ii(ll)});
  end
  fprintf(1,'---\n');
  ss = sum(apmap(:,ii),2);
  aa = find(ss);
  [ss ii] = sort(ss(aa));
  aa = aa(ii);
  for ll = length(aa):-1:1
    fprintf(1,'%d\t %s\n',ss(ll),authors{aa(ll)});
  end
  fprintf(1,'---\n');
  qq = struct(dpm.qq{kk});
  [ii jj mi]=find(qq.mi);
  [mi l] = sort(mi);
  for ll=length(mi):-1:max(1,length(mi)-10)
    fprintf(1,'%d\t %s\n',mi(ll),words{ii(l(ll))});
  end
  fprintf(1,'\n');
end
