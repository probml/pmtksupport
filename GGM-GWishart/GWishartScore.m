function [score, GWpost] = GWishartScore(X, G, delta0, S0, score_method)
%
% [score, GWpost] = GWishartScore(X, G, delta0, S0, Ptol, score_method)
%
% Description: 
%  This function implements a common interface to several different methods for scoring
%  Gaussian Graphical Models (GGMs) for use in stochastic local search or other applications.
%  The available methods include the BIC score, Laplace approximation to the marginal 
%  likelihood, diagonal Hessian Laplace approximation to the marginal likelihood, and 
%  a Monte Carlo approximation to the marginal likelihood as described in [1]. All of the 
%  methods except for the Monte Carlo approximation also estimate the MAP precision matrix P.
%  Estimation of the MAP precision matrix P is performed using Mark Schmidt's minFunc optimizer,
%  and is based on his sparse precision matrix estimation example. Detail are available at [4]. 
%  The Monte Carlo estimation code is based on code from Mike West's group at Duke. 
%  Details are available at [5].
%
% Inputs: 
%  X is n * d data matrix, presumed centered
%  G is d * d adjacency matrix
%  delta0 is prior strength; default 3
%  S0 is prior scatter matrix; default (delta0-2)*diag(cov(X,1))
%  score_method is one of {'mc','bic','laplace','diaglaplace'}
%    mc: Monte Carlo method of Atay-Kayis and Massam [2]. Sampling approximation to marginal likelihood
%    bic: loglik - dof/2 * log(N), dof = nedges + d
%    laplace: log p(D|G) approx -nd/2*log(2pi) + log Laplace(postZ) - log Laplace(priorZ)
%             This is the similar to  Lenkoski & Dobra [3], but they use Monte Carlo
%             for the priorZ term.
%    diaglaplace: same as laplace put only uses the diagonal elements of the Hessian
%
% Outputs:
%  score is the score of G given X and the prior as computed by score_method
%  GWpost contains additional information about the MAP estimate of P (computed by bic, laplace, and diaglaplace)
%
% -------------------------------------------------------------------------
% [1] Baback Moghaddam, Benjamin M. Marlin, Emtiyaz Khan, Kevin Murphy. Accelerating Bayesian 
%     Structural Inference for Non-Decomposable Gaussian Graphical Models. Proceedings of 
%     the Neural Information Processing Systems Conference (NIPS), 2009
% [2] Atay-Kayis & Massam. A Monte Carlo Method for Computing the Marginal likelihood in 
%     Non-decomposable Gaussian Graphical Models. Biometrika. 2005.
% [3] Lenkoski & Dobra, Bayesian structural learning & estimation in GGMs,
%     University of Washington, Dept. of Statistics, TR-545, 2009.
% [4] http://www.cs.ubc.ca/~schmidtm/Software/minFunc/minFunc.html 
% [5] http://www.stat.duke.edu/research/software/west/ggm.html
% -------------------------------------------------------------------------
%
% Revision History:
%  Baback Moghaddam 05/30/2009
%  Kevin Murphy 6/1/2009
%  Benjamin Marlin 6/21/2010

%Set default parameters
if nargin < 1, error('GWishartScore: Function GWishartScore requires and nxp data matrix as input');end
if nargin < 2, delta0 = 3; end
if nargin < 3, S0 = (delta0-2)*diag(diag(cov(X,1))); end
if nargin < 4, score_method='mc'; end

%Initialize the Data structure
Data.X    = X;
Data.XX   = X'*X;
[n p]     = size(Data.X);
noData.X  = zeros(0,p);  
noData.XX = 0;

%Initialize outputs
P      = [];
loglik = [];
GWpost = [];

%Compute MCMC estimator using mex code
if strcmpi(score_method,'mc')
  score   = ggm_gwish_mc_marglik(Data.X,Data.XX,G,S0,delta0,10,10000);
  return;
end

%Initialize GW prior structure
GWprior.d0 = delta0;
GWprior.S0 = S0;
GWprior.lognormconst = 0;
GWprior.lognormconstDiag = 0;

%If method is laplace, compute laplace approximation to denominator
if strcmpi(score_method,'laplace') | strcmpi(score_method,'diaglaplace')
  P0     = GWishartFit(noData, G, GWprior);
  GWtemp = compute_score(noData, G, P0, GWprior,score_method);
  GWprior.lognormconst = GWtemp.lognormconst;
  clear P0 GWtemp;
end
 
%Compute the map precision matrix P
P = GWishartFit(Data, G, GWprior);

%Compute the score 
[GWpost] = compute_score(Data, G, P, GWprior, score_method);
score    = GWpost.score;
loglik   = GWpost.loglik;
 
end


function [GWpost] = compute_score(Data, G, P, GWprior, score_method)
%
% Description:
%
% Finds the mode (MAP) estimate of GGM precision P given a Data struct and
% a conjugate G-Wishart prior W(d0,S0) struct GWprior (with fields: .d0 &.S0}
% and finds mode P of the G-Wishart posterior W(dn,Sn) with Sn=(dn-2)*inv(P)
% and returns the correponding (MAP) loglik and descired score.
% The score can be the BIC score, Laplace approximation to the marginal likelihood
% or a diagonal hessian approximation to the marginal likelihood. 
%
% Inputs:
%   Data.X  is an  n-by-p data matrix (presumed) centered
%   Data.XX is a  p-by-p scatter matrix (Data.X'*Data.X)
%   G is the GGM graph
%   P is the map estimate of the precision matrix under the GWishart prior as computed
%     by GWishartFit.
%   GWprior.d0  is  G-Wishart prior degree d0
%   GWprior.S0  is  G-Wishart p-by-p prior scatter matrix S0
%   GWprior.lognormconst  is  this G-Wishart's lognormconst (can be 0)
%   score_method is the score to compute {'bic','laplace','diaglaplace'}
%   Ptol is the tolerance for HTF covsel convergence (default = 1e-4)
%
% Outputs:
%   See code below for details of output structure GWposterior
% -------------------------------------------------------------------------

% Baback Moghaddam  05/30/2009
% Benjamin Marlin   06/21/2010

%Get data size
[n p] = size(Data.X);

%Set default tolerance
if nargin < 5
  Ptol = 1e-4; 
end

% get G-Wishart prior params
d0 = GWprior.d0;
S0 = GWprior.S0;

%Check that graph size matches number of data variables
if (size(G,1) ~= size(G,2)) | (size(G,1) ~= p)
  error('GWishartScore::Score: G must be p-by-p, with p = # cols in Data.X');
end

% check prior size violations
if (size(S0,1) ~= size(S0,2)) | (size(S0,1) ~= p)
  error('GWishartScore::Score: GWprior.S0 must be p-by-p, with p = # cols in Data.X \n');
end

% find W(dn,Sn) posterior mode P using HTFcovsel
dn = n + d0;
C = (Data.XX + S0) / (dn - 2);

% need logdetP and invP
[U,S,V] = svd(P);
invP = V*diag(1./diag(S))*U';
logdetP = sum(log(diag(S)));

% compute loglik
loglik   = -P(:)'*Data.XX(:)/2 + n*(logdetP - p*log(2*pi))/2;
numedge  = nnz(triu(G));
numdof   = numedge + p;
pcor     = cov2cor(P);

% the posterior Sn parameter
Sn = (dn-2)*invP;
logh = -(dn-2)*p/2 + (dn-2)*logdetP/2;  

% find full param set V
[Vi,Vj] = find(triu(G)+eye(p));
numdof = nnz(triu(G)+eye(p));
M1 = zeros(p);
M2 = zeros(p);

switch(score_method)
  %Compute BIC score
  case 'bic'
    if n > 0
      score = loglik - numdof*log(n)/2;
    else
      score = 0;  
    end
    GWpost.score   = score;

  %Diagonal hessian laplace approximation
  case 'diaglaplace'
    diagH = zeros(1, numdof);
    for e1 = 1:numdof
      e2 = e1;
      M1(:)=0; M1(:,[Vi(e1) Vj(e1)]) = invP(:,[Vj(e1) Vi(e1)]);
      M2(:)=0; M2(:,[Vi(e2) Vj(e2)]) = invP(:,[Vj(e2) Vi(e2)]);
      nz1 = [Vi(e1) Vj(e1)];
      nz2 = [Vi(e2) Vj(e2)];
   	  A = M1(nz2,nz1); B = M2(nz1,nz2);
	  tmp2 = A(1,:)*B(:,1) + A(2,:)*B(:,2);
      diagH(e1) = -(dn-2) * tmp2/2;
      %diagH(e1) = -(dn-2) * trace(M1(nz2,nz1)*M2(nz1,nz2))/2;
    end
    logdetHdiag          = sum(log(-diagH));
    lognormconst         = numdof*log(2*pi)/2 + logh - logdetHdiag/2;
    logmarglik           = lognormconst - GWprior.lognormconst - n*p*log(2*pi)/2;
    GWpost.lognormconst  = lognormconst;
    GWpost.score         = logmarglik;

  %Full laplace approximation
  case 'laplace'
    H = zeros(numdof);
    for e1 = 1:numdof
      for e2 = e1:numdof
          M1(:)=0; M1(:,[Vi(e1) Vj(e1)]) = invP(:,[Vj(e1) Vi(e1)]);
          M2(:)=0; M2(:,[Vi(e2) Vj(e2)]) = invP(:,[Vj(e2) Vi(e2)]);
          nz1 = [Vi(e1) Vj(e1)];
          nz2 = [Vi(e2) Vj(e2)];
          %tmp = trace(M1(nz2,nz1)*M2(nz1,nz2));
          A = M1(nz2,nz1); B = M2(nz1,nz2);
          tmp2 = A(1,:)*B(:,1) + A(2,:)*B(:,2);
          %assert(approxeq(tmp, tmp2));
          H(e1,e2) = -(dn-2) * tmp2/2;
          H(e2,e1) = H(e1,e2);
      end
    end
    logdetH = 2*sum(log(diag(chol(-H))));   % neg Hessian will be posdef
    lognormconst        = numdof*log(2*pi)/2 + logh - logdetH/2;
    logmarglik          = lognormconst - GWprior.lognormconst - n*p*log(2*pi)/2;
    GWpost.lognormconst = lognormconst;
    GWpost.score        = logmarglik;

end

% assign output struct
GWpost.dn      = d0 + n;
GWpost.P       = P;
GWpost.numedge = numedge;
GWpost.numdof  = numdof;
GWpost.logdetP = logdetP;
GWpost.loglik  = loglik;

end





