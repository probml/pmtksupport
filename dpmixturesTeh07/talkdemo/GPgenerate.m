% generate a sample function from a GP iteratively 
% (cheats by drawing whole function first)
ls = 2;
nn = 100;

xx = (0:nn)'/nn*10;

x2 = xx.^2;
x2 = x2(:,ones(1,nn+1));
vv = exp(-.5*(x2 + x2' - 2*xx*xx')/(ls^2));
[u s v] = svd(vv);

yy = u*(sqrt(diag(s)).*randn(nn+1,1));

clf;
axis([0 10 -2 2]);
set(gcf,'doublebuffer','on');

ii = randperm(nn+1);
aa = [0 10 floor(min(yy))-.5 ceil(max(yy))+.5];
for jj = [1:10 20:10:80]
  plot(xx(ii(1:jj)),yy(ii(1:jj)),'o','linewidth',3);
  title('A draw from a Gaussian process');
  axis(aa);
  fprintf(1,'\r%d out of [1:10 20:10:80]',jj);
  pause
  if jj<10,
    hold on;
    plot(xx(ii(jj+1)),0,'r+');
    hold off;
    pause
  end
end
hold on
plot(xx,yy,'k','linewidth',3);
hold off
axis(aa);
fprintf(1,'\n');
