function [f] = groupL1regularizer(w,lambda,groups)
nGroups = max(groups);
f = 0;
for g = 1:nGroups
    nCols = sqrt(sum(groups==g));
    f = f + lambda(g)*sum(svd(reshape(w(groups==g),nCols,nCols)));
end
end