function ll = logmarginal(qq);
% ll = logmarginal(qq)
% log marginal probability of data items in the component
% log p(x_1,...,x_n)

[ii jj mi] = find(qq.mi);
ll = qq.Z0 + gammaln(qq.aa*qq.dd) - gammaln(qq.aa*qq.dd+qq.mm) ...
           - length(mi)*gammaln(qq.aa) + full(sum(gammaln(qq.aa+mi)));
