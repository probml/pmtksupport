function dpm_demo2d_plot(dpm,cc);
% plots data and ellipses, assumes in 2D

ih = ishold;

for kk = 1:dpm.KK
  [mu sigma] = rand(dpm.qq{kk});
  plotellipse(mu,sigma,'color',cc(kk,:),'linewidth',2);
  hold on
  ii = find(dpm.zz==kk);
  xx = cat(2,dpm.xx{ii});
  plot(xx(1,:),xx(2,:),'.','color',cc(kk,:),'markersize',10);
end

if ~ih, hold off; end

