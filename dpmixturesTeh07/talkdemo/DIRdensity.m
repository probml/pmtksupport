% plots out density for Dirichlet distribution

aa = [5 2 2];
nn = 1001;
tiny = .0001;
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
lz = gammaln(sum(aa))-sum(gammaln(aa));
jj = logical(zeros(1,nn*nn));
for ii = 1:nn*nn
  if xx(ii) > tiny & yy(ii) >tiny & zz(ii) > tiny
    pp(ii) = -exp(log(xx(ii))*(aa(1)-1)+log(yy(ii))*(aa(2)-1)...
        +log(zz(ii))*(aa(3)-1)+lz);
  else
    pp(ii) = 0;
    jj(ii) = 1;
  end
end
%pp(jj) = min(
%pp = pp / max(pp) * .5;

colormap hot
%hh = contour(0:dd:1,0:dd:1,reshape(pp,[nn nn]),[5 5 5]);
%set(hh(1),'linestyle','none');
imagesc(reshape(pp,[nn nn]),[-10 0]);
hold on
plot([1.5 1.5 sqrt((nn-1)^2-(nn/2)^2)+.5 1.5],[.5 (nn-.5) (nn/2+.5) 0.5],'k');
hold off
axis([-1 sqrt((nn-1)^2-(nn/2)^2)+.5 -1 nn])
axis off
title(sprintf('Dir(%.1f,%.1f,%.1f)',aa(1),aa(2),aa(3)));
%axis([0 1 0 1]);
%axis([0 1 0 1 0 2.5]);
%set(gca,'CameraPosition', [-0.0464037 7.44273 8.0418]);
%set(gca,'CameraTarget', [0.5 0.5 0.25]);
%set(gca,'CameraUpVector', [0.1430892 -1.820793 0.40865]);
%set(gca,'CameraViewAngle', [4.6627]);
print('-depsc','-zbuffer',...
        ['dirc-' num2str(aa(1)) '-' num2str(aa(2)) '-' num2str(aa(3)) '.eps']);




