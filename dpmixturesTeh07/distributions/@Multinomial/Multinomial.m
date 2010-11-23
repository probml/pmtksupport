function qq = Multinomial(hh);
% qq = Multinomial(hh)
% Creates a multinomial component containing no data items,
% and symmetric Dirichlet prior given by hh.
%
% MODEL:
% hh.dd         = (1x1) dimensionality
% hh.aa         = (1x1) parameter for symmetric Dirichlet
% qq.pi         = (dx1) multinomial probabilities
% xx            = (dx1) sparse counts
% mm            = (1x1) total counts in xx
%
% qq.pi         ~ Dirichlet(hh.aa/hh.dd,...,hh.aa/hh.dd)
% xx | mm,hh,qq ~ Multinomial(mm,qq.pi)

qq.dd = hh.dd;
qq.aa = hh.aa/hh.dd;
qq.mi = sparse(hh.dd,1);
qq.mm = 0;
qq.nn = 0;
qq.Z0 = 0;
qq    = class(qq,'Multinomial');
