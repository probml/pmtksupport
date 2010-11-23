clear
dd = 1;
s0 = 3;
ss = 1;

hh.dd = dd;
hh.ss = s0^2/ss^2;
hh.vv = 5;
hh.VV = ss^2*eye(dd);
hh.uu = zeros(dd,1);

NN = 10;
xx = [-2+.5*randn(1,5) 2+.5*randn(1,5)];
yy = -10:.01:10;

qq = Gaussian(hh)
clf
hold on
for iter = 1:10
  [mu sigma] = rand(qq);
  plot(yy,1/sqrt(2*pi*sigma)*exp(-.5/sigma*(yy-mu).^2));
end
pause

for ii = 1:10
  qq = additem(qq,xx(ii));
  clf
  plot(xx(1:ii),zeros(1,ii),'kx','linewidth',5);
  hold on
  for iter = 1:10
    [mu sigma] = rand(qq);
    plot(yy,1/sqrt(2*pi*sigma)*exp(-.5/sigma*(yy-mu).^2));
  end
  pause
end


