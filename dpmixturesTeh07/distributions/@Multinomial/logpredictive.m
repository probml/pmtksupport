function ll = logpredictive(qq,xx);
% ll = logpredictive(qq,xx)
% log predictive probability of xx given other data items in the component
% log p(xx|x_1,...,x_n)

if issparse(xx)
  [ii jj mi] = find(xx);
  mm = sum(mi);
  ll = gammaln(mm+1) - sum(gammaln(mi+1)) ...
       + gammaln(qq.aa*qq.dd+qq.mm) - gammaln(qq.aa*qq.dd+qq.mm+mm) ...
       + full(sum(gammaln(qq.aa+qq.mi+xx) - gammaln(qq.aa+qq.mi)));
elseif isscalar(xx)
  ll = log((qq.aa+qq.mi(xx))/(qq.aa*qq.dd+qq.mm));
else
  error('data item xx type unknown.');
end
