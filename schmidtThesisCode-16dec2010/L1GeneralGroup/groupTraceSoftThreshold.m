function [w] = groupTraceSoftThreshold(w,alpha,lambda,groups)
    nGroups = max(groups);
    for g = 1:nGroups
        gNdx = groups==g;
        nCols = sqrt(sum(gNdx));
        W = reshape(w(groups==g),nCols,nCols);
        [U,S,V] = svd(W);
        S = sign(S).*max(0,S-lambda(g)*alpha);
        W = U*S*V';
        w(groups==g) = W(:);
    end
end