function ll = logmarginal(qq);
% ll = logmarginal(qq)
% log marginal probability data items in the component
% log p(x_1,...,x_n)

ll = ZZ(qq.dd,qq.nn,qq.rr,qq.vv,qq.CC,qq.XX) - qq.Z0;

