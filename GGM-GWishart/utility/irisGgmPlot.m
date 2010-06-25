% Plot structure of 4 node graph on iris virginica data
% Same as Atay-Kayis & Massam 2005 fig 3
function plotIrisGgm(G,showlabels)

% SL SW PL PW
tx = [-0.8 1.1 -0.8 1.1];
ty = [1.1 1.1 -0.2 -0.2];
x  = [0 1 0 1];
y  = [1 1 0 0];
names = {'SL', 'SW', 'PL', 'PW'};
if(showlabels)
  for i=1:4
    text(tx(i),ty(i),sprintf('%s',names{i}));
  end
end
hold on
for i=1:4
  for j=i+1:4
    if G(i,j)
      plot([x(i) x(j)], [y(i) y(j)],'.b-','MarkerFaceColor','b'); hold on;
    end
  end
end

axis([-0.5 1.1 -0.2 1.1])
axis off
end