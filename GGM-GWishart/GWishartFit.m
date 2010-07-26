function P = GWishartFit(Data, G, GWprior)
%
% P = GWishartFit(Data, G, GWprior)
%
% Description:
% Finds the mode (MAP) estimate of GGM precision matrix P given a Data struct and
% a conjugate G-Wishart prior W(d0,S0) struct GWprior (with fields: .d0 &.S0)
% The mode finding algorithm simply optimizes the GGM log likelihood with respect 
% to the free elements of P. This is a convex optimization problem. We solve it 
% using Mark Schmidt's minFunc optimizer, which must be installed to use this code. 
%
% Inputs:
%   Data.X  is an  n-by-p data matrix (presumed) centered
%   Data.XX is a  p-by-p scatter matrix (Data.X'*Data.X)
%   GWprior.d0  is  G-Wishart prior degree d0
%   GWprior.S0  is  G-Wishart p-by-p prior scatter matrix S0
%
% Outputs:
%   P is the MAP precision matrix under the GWishart prior
% -------------------------------------------------------------------------
%
% Revision History:
%   Baback Moghaddam  05/30/2009
%   Benjamin Marlin   06/21/2010

%Get data size
[n p] = size(Data.X);

% get G-Wishart prior params
d0 = GWprior.d0;
S0 = GWprior.S0;

%Check that graph size matches number of data variables
if (size(G,1) ~= size(G,2)) || (size(G,1) ~= p)
  error('GWishartScore::Fit: G must be p-by-p, with p = # cols in Data.X');
end

% check prior size violations
if (size(S0,1) ~= size(S0,2)) || (size(S0,1) ~= p)
  error('GWishartScore::Fit: GWprior.S0 must be p-by-p, with p = # cols in Data.X \n');
end

% compute posterior scatter matrix
dn = n + d0;
C = (Data.XX + S0) / (dn - 2);

%Find non-zero elements of upper triangle of G
%make sure diagonal is non-zero
nonZero = triu((eye(p)+G)>0);
nonZero = nonZero(:);

%Set minfunc options
options.TolX            = 1e-16;           
options.TolFun          = 1e-16;
options.Method          = 'lbfgs';
options.MaxFunEvals     = 5000;
options.MaxIter         = 5000;    
options.DerivativeCheck = 'off';
options.Display         = 'off';

%Define the objective function, run the optimizer and get the results
P0         = eye(p);
lambda     = 0;
funObj     = @(x)sparsePrecisionObjLambda(x,p,nonZero,C,lambda);
[Ptmp,f]   = minFunc(funObj,P0(nonZero),options);
P          = zeros(p);
P(nonZero) = Ptmp;
P          = P + triu(P,1)';

end

function [f,g] = sparsePrecisionObjLambda(x,p,nonZero,C,lambda)

%Description: This function compute the objective and gradient for
%  the MLE of the precision matrix given emperical covariance matrix C
%  and list of non-zero upper triangle precision matrix entries given by
%  nonZero. lambda is an l2 regularizer on the precision matrix entries. 
%  Based on sparse GGM estimation code by Mark Schmidt.
% Revision History:
%   Benjamin Marlin   06/21/2010

  X = zeros(p);  %initialize X to zeros
  X(nonZero) = x;    %fill the diagonal and upper triangle
  X = X + triu(X,1)';%fill the lower triangle

  %compute cholesky factorization to check if current precision matrix
  %is positive definite
  [R,posdefind] = chol(X); 

  if(posdefind == 0)
      % Matrix is in positive-definite cone
      % Fast Way to compute -logdet(X) + tr(X*C)
      f = -2*sum(log(diag(R))) + sum(sum(C.*X)) + (lambda/2)*sum(X(:).^2);
      if(nargout>1)
        g = -inv(X) + C + lambda*X; %compute gradient
        g = g + tril(g,-1)'; %add contribution from lower triangle to upper triangle
        g = g(nonZero); %retain diagonal and upper triangle only
      end
  else
      % Matrix not in positive-definite cone, set f to Inf
      % to force minFunc to backtrack
      f = inf;
      g = 0;
  end

end
