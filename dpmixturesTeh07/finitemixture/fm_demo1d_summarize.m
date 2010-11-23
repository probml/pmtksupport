function fm_demo1d_summarize(fm,yy,pp);
% summarize posterior densities, assumes in 1D

clf
hold on

% plot data items
xx = cat(2,fm.xx{:});
plot(xx,zeros(size(xx)),'kx','markersize',3,'linewidth',3);

nump = size(pp,1);
meanp = mean(pp,1);
sortp = sort(pp,1);
q5p = mean(sortp(round(nump/100*4):round(nump/100*6),:),1); % 5th quantile
q50p = mean(sortp(round(nump/100*49):round(nump/100*51),:),1); % 50th quantile
q95p = mean(sortp(round(nump/100*94):round(nump/100*96),:),1); % 95th quantile

% plot 5-95 quantiles of densities
patch([yy fliplr(yy)],[q5p fliplr(q95p)],[.9 .9 .9],'edgecolor',[.9 .9 .9]);


% plot mean
plot(yy,meanp,'r','linewidth',3);

% plot median
plot(yy,q50p,'b-.','linewidth',2);

% plot 5 samples
rr = randperm(nump);
plot(yy,pp(rr(1:5),:));

legend('5-95 quantile','mean','median');
hold off
