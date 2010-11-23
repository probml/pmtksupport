function qq = Gaussian(hh);
% qq = Gaussian(hh)
% Creates a gaussian component containing no data items, 
% and gaussian-wishart prior given by hh.
%
% MODEL:
% hh.dd       = (1x1) dimensionality
% hh.ss       = (1x1) relative variance of mm versus data (cluster separability)
% hh.VV       = (dxd) mean covariance matrix of clusters.
% hh.vv       = (1x1) degrees of freedom of inverse Wishart covariance prior.
% hh.uu       = (dx1) prior mean vector
% qq.RR       = (dxd) precision matrix
% qq.mm       = (dx1) mean vector
% xx          = (dx1) data vector
%
% qq.RR       ~ Wishart(hh.dd, hh.vv, hh.SS)
% qq.mm       ~ Normal(hh.dd, hh.uu, hh.rr*qq.RR)
%
% xx | hh, qq ~ Normal(hh.dd, qq.mm, qq.RR) iid
%

hh.rr = 1/hh.ss;
hh.SS = hh.VV*(hh.vv); %-hh.dd-1);

qq.dd = hh.dd;
qq.nn = 0;                  % number of items.
qq.rr = hh.rr;
qq.vv = hh.vv;
qq.CC = chol(hh.SS + hh.rr*hh.uu*hh.uu');
qq.XX = hh.rr*hh.uu;
qq.Z0 = ZZ(hh.dd,qq.nn,qq.rr,qq.vv,qq.CC,qq.XX);

qq = class(qq,'Gaussian');
