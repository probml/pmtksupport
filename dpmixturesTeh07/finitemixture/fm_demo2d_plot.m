function fm_demo2d_plot(fm,cc);
% plots data and ellipses, assumes in 2D

ih = ishold;

for kk = 1:fm.KK
  [mu sigma] = rand(fm.qq{kk});
  plotellipse(mu,sigma,'color',cc(kk,:));
  hold on
  ii = find(fm.zz==kk);
  if length(ii) > 0
    xx = cat(2,fm.xx{ii});
    plot(xx(1,:),xx(2,:),'.','color',cc(kk,:));
  end
end

if ~ih, hold off; end

