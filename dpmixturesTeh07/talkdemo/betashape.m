aa = [.7 1 2 3.5 5];

xx = 0:.01:1;

yy = zeros(length(aa),length(xx));
for ii = 1:length(aa)
  yy(ii,:) = betapdf(xx,1,aa(ii));
end
legend('\alpha=.7','\alpha=1','\alpha=2','\alpha=3.5','\alpha=5');

hh = plot(xx,yy');
set(hh,'linewidth',2);
title('pdf of Beta(1,\alpha)');

print -depsc betashape
