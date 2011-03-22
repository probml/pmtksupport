function [loglik,g] = studentLoss(X,y,w,lambda,dof)
[nInstances,nVars] = size(X);

if lambda > 0 && dof > 0
    logZ = gammaln(dof/2+1/2) - gammaln(dof/2) + (1/2)*log(lambda) - (1/2)*log(pi*dof);
    err = X*w-y;
    logLiks = 1 + (lambda/dof)*err.^2;
    sumLL = sum(log(logLiks));
    loglik = -nInstances*logZ - (-dof/2 - 1/2)*sumLL;
    
    if nargout > 1
        g = zeros(nVars+2,1);
        g(1:nVars) = -2*(lambda/dof)*(-dof/2 - 1/2)*X'*(err./logLiks);
        g(nVars+1,1) = -(nInstances/2)*(1/lambda) - (-dof/2 - 1/2)*(1/dof)*sum((err.^2)./logLiks);
        g(end,1) = -(nInstances/2)*psi(dof/2+1/2) + (nInstances/2)*psi(dof/2) + nInstances/(2*dof);
        g(end,1) = g(end,1) + (-1/2 - dof/2)*(lambda/dof^2)*sum((err.^2)./logLiks);
        g(end,1) = g(end,1) + (1/2)*sumLL; 
    end
else
    loglik = inf;
    g = zeros(nVars+2,1);
end
end