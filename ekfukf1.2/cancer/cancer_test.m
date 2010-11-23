load /proj/finnwell/spatial/data/cancer2007/testdata/lungcancer_female_10years.mat
%name = sprintf('c%d_lungcancer_FIC_3exp_.mat',run);


[n, nin] = size(xx);

% $$$    gpcf1 = gpcf_st_exp('init',2,1, 'lengthScale', [20 26], 'magnSigma2', 0.04);
% $$$    gpcf1.p.lengthScale = t_p({1 4});
% $$$    gpcf1.p.magnSigma2 = t_p({0.3 4});
% $$$ 
% $$$    gpcf2 = gpcf_t_exp('init',3, 'lengthScale', [100], 'magnSigma2', 0.11);
% $$$    gpcf2.p.lengthScale = t_p({1 4});
% $$$    gpcf2.p.magnSigma2 = t_p({0.3 4});
% $$$ 
% $$$    gpcf3 = gpcf_s_exp('init',3, 'lengthScale', [90], 'magnSigma2', 0.11);
% $$$    gpcf3.p.lengthScale = t_p({1 4});
% $$$    gpcf3.p.magnSigma2 = t_p({0.3 4});

load '/proj/finnwell2/jmjharti/ajot/cancer/lungcancer_FIC/c51_lungcancer_FIC_3exp_.mat'
rt = thin(rgp, 1000, 80);
   
nsamp = length(rgp.hmcrejects);
subplot(7,1,1)
plot(1:nsamp,rgp.cf{1}.lengthScale(:,1));
subplot(7,1,2)
plot(1:nsamp,rgp.cf{1}.lengthScale(:,2));
subplot(7,1,3)
plot(1:nsamp,rgp.cf{1}.magnSigma2);
subplot(7,1,4)
plot(1:nsamp,rgp.cf{2}.lengthScale(:,1));
subplot(7,1,5)
plot(1:nsamp,rgp.cf{2}.magnSigma2);
subplot(7,1,6)
plot(1:nsamp,rgp.cf{3}.lengthScale(:,1));
subplot(7,1,7)
plot(1:nsamp,rgp.cf{3}.magnSigma2);
drawnow

xx2 = xx(1:431,1:2);
rgp2 = 