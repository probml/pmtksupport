function l = ggm_gwish_mc_marglik(X,scatter,graph,S0,deltaGW,ntrials,nsamples)
%
% l = ggm_gwish_mc_marglik(scatter,graph,S0,delta0)
%
% Inputs:
%   X       - the data set with data cases as rows
%   scatter - the empirical scatter matrix. Will be computed if
%             an empty matrix is passed in.
%   graph   - the adjacency matrix for the GGM graph
%   S0      - the scale parameter for the gWishart distribution
%   deltaGW - the degrees of freedom parameter from the gWishart distribution
%   ntrials - maximum number of sampling trials to perform
%   nsamples- number of samples to draw per trial
%
% Outputs:
%   l       - the marginal likelihood of X under G.
%
%Details: This function is a wrapper for the mex file
%mex_ggm_hiw_mc_marglik, which is an implementation of
%the Atay-Kayis/Massam MC method for computing the
%marginal likelihood of a data set under a Gaussian likelihood
%and hyper inverse Wishart prior. The mex file is an adaptation and
%extension of stand-alone C++ code made available by Mike West's group 
%at Duke. The original C++ code is available here:
%http://www.stat.duke.edu/research/software/west/ggm/MH-u.tar.gz
%Note that we convert from the GW parametrization to the
%HIW parametrization using deltaHIW = deltaGW-p+1 where p
%is the data dimension.
%
%Author: Benjamin Marlin
%Date:   Sept 18,2009
%

%Get number of data cases and data dimensions
[n,p] = size(X);

%Make scatter matrix of empty passed in
if(isempty(scatter))
  scatter = cov(X,1)*n;
end

%Convert the graph to an edge list
%Make sure the diagonal is non-zero
[g(:,1),g(:,2),foo] = find(graph+eye(p));
g                   = int32(g);

%Check S0 argument. Only spherical S0 are allowed (ie: S0 = s*I)
if( ~all(all(S0 == S0(1)*eye(p))))
  error('ggm_gwish_mc_marglik: Currently only spherical S0 parameters are supported (ie: S0=s*I)');
end

%Call mex version of West code for Atay-Kayis/Massam MC method 
%marginal likelihood under the HIW prior. Final argument in the
%call is the initial state of the random number generator used in the
%mex code. For more info view WEST-README and source code for mex file.
l = mex_ggm_hiw_mc_marglik(scatter,g,n,S0(1),deltaGW-p+1,int32(ntrials),int32(nsamples),uint32((2^31)*rand));
if(isnan(l));l=-1e100;end;