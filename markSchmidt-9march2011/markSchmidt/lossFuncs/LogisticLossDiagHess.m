function [nll,g,D] = LogisticDiagHessLoss(w,X,y)
% w(feature,1)
% X(instance,feature)
% y(instance,1)

[n,p] = size(X);

Xw = X*w;
yXw = y.*Xw;

nll = sum(mylogsumexp([zeros(n,1) -yXw]));

if nargout > 1
    if nargout > 2
        sig = 1./(1+exp(-yXw));
        g = -(X.'*(y.*(1-sig)));
    else
        %g = -X.'*(y./(1+exp(yXw)));
        g = -(X.'*(y./(1+exp(yXw))));
    end
end

if nargout > 2
    % Compute diagonals of Hessian
    sig = sig.*(1-sig);
%     D = zeros(p,1);
%     for i = 1:p
%         D(i,1) = (sig.*X(:,i))'*X(:,i);
%     end
    D = sum(diag(sparse(sig))*X.^2,1)';
    D(D==0)=inf; % Hack
end