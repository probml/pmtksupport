function qq = delitem(qq,xx);
% qq = delitem(qq,xx)
% deletes data item xx into component qq.
% xx can either be a sparse dx1 vector of counts,
% or a scalar indicating value of a single draw.

if issparse(xx)
  [ii jj mi] = find(xx);
  mm    = sum(mi);
  qq.nn = qq.nn - 1;
  qq.mi = qq.mi - xx;
  qq.mm = qq.mm - mm;
  qq.Z0 = qq.Z0 - gammaln(mm+1) + sum(gammaln(mi+1));
elseif isscalar(xx)
  qq.nn     = qq.nn - 1;
  qq.mi(xx) = qq.mi(xx) - 1;
  qq.mm     = qq.mm - 1;
else
  error('data item xx type unknown.');
end

