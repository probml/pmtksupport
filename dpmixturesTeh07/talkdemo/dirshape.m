aa = [.9 .9 .9];
nn = 101;
dd = 1/(nn-1);

MM = [cos(pi/3) 1 0; sin(pi/3) 0 0; 1 1 1];
IM = inv(MM);

x1 = reshape((0:dd:1)'   * ones(1,nn), 1, nn*nn);
y1 = reshape(ones(1,nn)' * (0:dd:1)  , 1, nn*nn);
z1 = reshape(ones(1,nn)' * ones(1,nn), 1, nn*nn);

XX = IM * [x1; y1; z1];
xx = XX(1,:);
yy = XX(2,:);
zz = XX(3,:);

pp = zeros(1,nn*nn);
for ii = 1:nn*nn
  xyz = [xx(ii) yy(ii) zz(ii)];
  jj = find(xyz>0);
  pp(ii) = exp(gammaln(sum(aa(jj)))-sum(gammaln(aa(jj)))+(aa(jj)-1)*xyz(jj)');
end
pp = pp / max(pp) * .5;

colormap default
mm = contour(0:dd:1,0:dd:1,reshape(pp,[nn nn]),20);
hold on
plot(MM(2,[1 2 3 1]),MM(1,[1 2 3 1]));
hold off
axis off
set(gca,'CameraPosition', [-0.0464037 7.44273 8.0418]);
set(gca,'CameraPositionMode', 'manual');
set(gca,'CameraTarget', [0.5 0.5 0.25]);
set(gca,'CameraTargetMode', 'auto');
set(gca,'CameraUpVector', [0.1430892 -1.820793 0.40865]);
set(gca,'CameraUpVectorMode', 'manual');
set(gca,'CameraViewAngle', [4.6627]);
set(gca,'CameraViewAngleMode', 'auto');
axis([0 1 0 1 0 .5])
%print -depsc dircontour

colormap gray
hh = surfl(0:dd:1,0:dd:1,reshape(pp,[nn nn]),[5 5 5]);
set(hh(1),'linestyle','none');
hold on
plot(MM(2,[1 2 3 1]),MM(1,[1 2 3 1]));
hold off
axis off
set(gca,'CameraPosition', [-0.0464037 7.44273 8.0418]);
set(gca,'CameraPositionMode', 'manual');
set(gca,'CameraTarget', [0.5 0.5 0.25]);
set(gca,'CameraTargetMode', 'auto');
set(gca,'CameraUpVector', [0.1430892 -1.820793 0.40865]);
set(gca,'CameraUpVectorMode', 'manual');
set(gca,'CameraViewAngle', [4.6627]);
set(gca,'CameraViewAngleMode', 'auto');
%print -depsc -zbuffer dirhill



