% Description: This example compares scoring functions that approximate the 
% marginal likelihood of a GGM. The comparison includes the BIC score,
% Laplace approximation (related to  Lenkoski & Dobra [3]), diagonal Hessian 
% Laplace approximation (discussed in Moghaddam et al. [1]), and the 
% Monte Carlo approximation of Atay-Kayis & Massam [2]. Implementation details
% and further sources can be found in GWishartScore.m and GWishartFit.m  
% For more information on these marginal likelihood approximations, see 
% Moghaddam et al. [1].
%
% This example computes the marginal likelihood approximation given by each of 
% the above methods for all four-node GGMs using the Iris virginica data. The 
% results can be compared to Figure 3 of Atay-Kayis & Massam [2, p 333]. 
%
% The BIC, Laplace and diagonal Hessian Laplace methods require Mark Schmidt's 
% free minFunc optimization package. minFunc can be downloaded from here:
% http://www.cs.ubc.ca/~schmidtm/Software/minFunc.html 
%
% -------------------------------------------------------------------------
% [1] Baback Moghaddam, Benjamin M. Marlin, Emtiyaz Khan, Kevin Murphy. Accelerating Bayesian 
%     Structural Inference for Non-Decomposable Gaussian Graphical Models. Proceedings of 
%     the Neural Information Processing Systems Conference (NIPS), 2009.
% [2] Atay-Kayis & Massam. A Monte Carlo Method for Computing the Marginal likelihood in 
%     Non-decomposable Gaussian Graphical Models. Biometrika. 2005.
% [3] Lenkoski & Dobra, Bayesian structural learning & estimation in GGMs,
%     University of Washington, Dept. of Statistics, TR-545, 2009.
% -------------------------------------------------------------------------
%
% Revision History:
%  Kevin Murphy 6/1/2009
%  Benjamin Marlin 6/21/2010

  %Set and check paths, packages, and compiled files
  addpath('utility');
  addpath('west-mc');
  if(isempty(which('minFunc')))
    error('gWishartExample: Please set your path to minFunc. See help for instructions on downloading the free minFunc optimization package.') 
  end
  if(exist('mex_ggm_hiw_mc_marglik')~=3)
    cd west-mc
    make_mex_ggm_hiw_mc_marglik
    cd ..
  end

  %Configuration options
  score_print_names = {'MC','BIC','Diag Laplace','Laplace'};
  score_methods     = {'mc','bic','diaglaplace','laplace'};
  nM                = length(score_methods);
  
  %Load iris data matrix
  load('data/irisData.mat');
  c      = find(strcmp(classnames, 'virginica'));
  ndx    = find(y==c);
  X      = X(ndx,:);
  X      = bsxfun(@plus,X,-mean(X,1));
  [n,p]  = size(X);

  %Make all undirected graphs
  G          = mk_all_ugs(p);
  num_graphs = length(G);   

  %Set GW priors
  d0  =  3+p-1;  % Roverato'02: convert from HIW to G-Wishart (delta + |V| - 1 = 3 + 4 - 1  )
  S0  =  eye(p); % same as Roverato'02

  %Compute scores for all graphs 
  score = zeros(num_graphs,nM);
  for g=1:num_graphs
    fprintf('Computing score for graph %d/%d\n',g,num_graphs); 
    for i=1:nM
      score(g,i) = GWishartScore(X, G{g}, d0, S0, score_methods{i}); 
    end
  end

  %Exponentiate and normalize scores to get graph posteriors
  pG = zeros(num_graphs,nM);
  for i=1:nM
    pG(:,i) = exp(normalizeLogspace(score(:,i)'))';
  end

  %Plot top graphs
  figure; 
  c=1; 
  for i=1:nM
    [foo,order] = sort(pG(:,i), 'descend');
    for k=1:8
      subplot(nM,8,c)
      irisGgmPlot(G{order(k)},c==1);
      title(sprintf('   %5.3f', pG(order(k),i)));
      if(k==1);
         text(-1.2,0.5,score_print_names{i},'Rotation',90,'HorizontalAlignment','center'); 
       end; 
      c=c+1; 
    end
  end
  set(gcf,'name','Top Graphs with Posterior Probabilities');

  
  %Compute and plot KL to MC solution
  KL = zeros(nM,1);
  for i=2:nM
    KL(i) = sum(pG(:,1).*(-log2(eps+pG(:,i)) + log2(eps+pG(:,1))),1);
  end
  figure;set(gcf,'name','KL Divergence to MC Posterior');
  bar(KL(2:end)); set(gca,'XTickLabel',score_print_names(2:end));
  ylabel('KL Divergence (bits)');
  xlabel('Scoring Method');
  title('KL Divergence from MC Posterior');
  grid on;
  