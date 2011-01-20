% THE GP TOOLS (/in the GP/ folder):
% 
%  Gaussian process utilities:
%   GPCOV     Evaluate covariance matrix between two input vectors. 
%   GP_TRCOV  Evaluate training covariance matrix (gp_cov + noise covariance). 
%   GP_TRVAR  Evaluate training variance vector. 
%   GP_PAK    Combine GP hyper-parameters into one vector.
%   GP_UNPAK  Set GP hyper-parameters from vector to structure
%   GP_RND    Random draws from the postrior Gaussian process
%   GP_INIT   Create a Gaussian Process data structure. 
%
%  Covariance functions:
%   GPCF_CONSTANT      Create a constant covariance function 
%   GPCF_EXP           Create a squared exponential covariance function
%   GPCF_LINEAR        Create a linear covariance function
%   GPCF_MATERN32      Create a squared exponential covariance function
%   GPCF_MATERN52      Create a squared exponential covariance function
%   GPCF_NEURALNETWORK Create a squared exponential covariance function
%   GPCF_NOISE         Create a noise covariance function 
%   GPCF_NOISET        Create a scale mixture noise covariance function (~Student-t) 
%   GPCF_PERIODIC      Create a periodic covariance function
%   GPCF_PPCS0         Create a piece wise polynomial (q=0) covariance function 
%   GPCF_PPCS1         Create a piece wise polynomial (q=1) covariance function 
%   GPCF_PPCS2         Create a piece wise polynomial (q=2) covariance function 
%   GPCF_PPCS3         Create a piece wise polynomial (q=3) covariance function 
%   GPCF_PROD          Create a product form covariance function 
%   GPCF_RQ            Create an rational quadratic covariance function 
%   GPCF_SEXP          Create a squared exponential covariance function
%
%  Likelihood functions:
%   LIKELIH_BINOMIAL   Create a binomial likelihood structure 
%   LIKELIH_LOGIT      Create a Logit likelihood structure 
%   LIKELIH_NEGBIN     Create a Negbin likelihood structure 
%   LIKELIH_POISSON    Create a Poisson likelihood structure 
%   LIKELIHOOD_PROBIT  Create a Probit likelihood structure 
%   LIKELIH_T          Create a Student-t likelihood structure 
%
% Inference utilities:
%   GP_E          Evaluate energy function (un-normalized marginal log posterior) 
%                 in case of Gaussian observation model
%   GP_G          Evaluate gradient of energy (GP_E) for Gaussian Process
%   GP_PRED       Make predictions with Gaussian process 
%   GPEP_E        Conduct Expectation propagation and return marginal 
%                 log posterior estimate
%   GPEP_G        Evaluate gradient of EP's marginal log posterior estimate (GPEP_E)
%   EP_PRED       Predictions with Gaussian Process EP approximation
%   GPLA_E        Construct Laplace approximation and return marginal 
%                 log posterior estimate
%   GPLA_G        Evaluate gradient of Laplace approximation's marginal 
%                 log posterior estimate (GPLA_E)
%   LA_PRED       Predictions with Gaussian Process Laplace approximation
%   GP_MC         Markov chain sampling for Gaussian process models
%   MC_PRED       Predictions with Gaussian Process MCMC approximation.
%   GP_IA         Integration approximation with grid, Monte Carlo or CCD integration
%   IA_PRED       Prediction with Gaussian Process GP_IA solution.
%    LGCP         Log Gaussian Cox Process intensity estimate for 1D and 2D data
%
%  Model assesment and comparison:
%   GP_DIC        The DIC statistics and efective number of parameters in a GP model
%   GP_KFCV       K-fold cross validation for a GP model
%   GP_LOOE       Evaluate the leave one out predictive density in case of
%                 Gaussian observation model
%   GP_LOOE       Evaluate the gradient of the leave one out predictive 
%                 density (GP_LOOE) in case of Gaussian observation model 
%   GP_PEFF       The efective number of parameters in GP model with focus 
%                 on latent variables.
%
%  Metrics:
%   METRIC_DISTANCEMATRIX  An Euclidean distance for Gaussian process models. 
%   METRIC_EUCLIDEAN       An Euclidean distance for Gaussian process models.
%   METRIC_IBS_GXE         An Euclidean distance for Gaussian process models.
%  
%  Misc:
%    LDLROWMODIFY  Function to modify the sparse cholesky factorization 
%                  L*D*L' = C, when a row and column k of C have changed 
%    LDLROWUPDATE  Multiple-rank update or downdate of a sparse LDL' factorization.
%    SPINV         Evaluate the sparsified inverse matrix
%    SCALED_HMC    A scaled hybric Monte Carlo samping for latent values
%    SCALED_MH     A scaled Metropolis Hastings samping for latent values
%    TRCOV         Evaluate training covariance matrix for covariance function
%                  This is a mex-function that is called from gpcf_*_trcov
%                  functions.
%    GP_INSTALL    Matlab function to compile all the c-files to mex in the 
%                  GPstuff/gp folder.
%
%  Demonstration programs:
%   DEMO_BINOMIAL          Demonstration of Gaussian process model with binomial
%                          likelihood
%   DEMO_BINOMIAL2         Demonstration for modeling age-period-cohort data
%                          by a binomial model combined with GP prior.
%   DEMO_CLAASIFIC         Classification problem demonstration for 2 classes 
%   DEMO_COMPARESPARSEGP   Regression demo comparing different sparse
%                          approximations
%   DEMO_LGCP              Demonstration for a log Gaussian Cox process
%                          with inference via EP or Laplace approximation
%   DEMO_MODELASSESMENT1   Demonstration for model assesment with DIC, number 
%                          of effective parameters and ten-fold cross validation
%   DEMO_MODELASSESMENT2   Demonstration for model assesment when the observation 
%                          model is non-Gaussian
%   DEMO_INFNEURALNETWORK  Demonstration of Gaussian process with a neural
%                          network covariance function
%   DEMO_PERIODICCOV       Regression problem demonstration for periodic data
%   DEMO_PPCSCOV           Regression problem demonstration for 2-input 
%                          function with Gaussian process using CS covariance
%   DEMO_REGRESSION1       Regression problem demonstration for 2-input 
%                          function with Gaussian process
%   DEMO_REGRESSION2       Regression problem demonstration with additive model
%   DEMO_REGRESSION_ADDITIVE Regression demonstration with additive Gaussian
%                          process using linear, squared exponential and
%                          neural network covariance fucntions 
%   DEMO_ROBUSTREGRESSION  A regression demo with Student-t distribution as a 
%                          residual model.
%   DEMO_SPARSEREGRESSION  Regression problem demonstration for 2-input 
%                          function with sparse Gaussian processes
%   DEMO_SPATIAL1          Demonstration for a disease mapping problem
%                          with Gaussian process prior and Poisson likelihood
%   DEMO_SPATIAL2          Demonstration for a disease mapping problem with 
%                          Gaussian process prior and negative binomial 
%                          observation model
% THE DIAGNOSTIC TOOLS (/in the diag/ folder):
%
% Covergence diagnostics
%   PSRF     - Potential Scale Reduction Factor
%   CPSRF    - Cumulative Potential Scale Reduction Factor
%   MPSRF    - Multivariate Potential Scale Reduction Factor
%   CMPSRF   - Cumulative Multivariate Potential Scale Reduction Factor
%   IPSRF    - Interval-based Potential Scale Reduction Factor
%   CIPSRF   - Cumulative Interval-based Potential Scale Reduction Factor
%   KSSTAT   - Kolmogorov-Smirnov goodness-of-fit hypothesis test
%   HAIR     - Brooks' hairiness convergence diagnostic
%   CUSUM    - Yu-Mykland convergence diagnostic for MCMC
%   SCORE    - Calculate score-function convergence diagnostic
%   GBINIT   - Initial iterations for Gibbs iteration diagnostic
%   GBITER   - Estimate number of additional Gibbs iterations
%
% Time series analysis
%   ACORR      - Estimate autocorrelation function of time series
%   ACORRTIME  - Estimate autocorrelation evolution of time series (simple)
%   GEYER_ICSE - Compute autocorrelation time tau using Geyer's
%                initial convex sequence estimator
%                (requires Optimization toolbox) 
%   GEYER_IMSE - Compute autocorrelation time tau using Geyer's
%                initial monotone sequence estimator
%
% Kernel density estimation etc.:
%   KERNEL1  - 1D Kernel density estimation of data
%   KERNELS  - Kernel density estimation of independent components of data
%   KERNELP  - 1D Kernel density estimation, with automatic kernel width
%   NDHIST   - Normalized histogram of N-dimensional data
%   HPDI     - Estimates the Bayesian HPD intervals
%
% Manipulation of MCMC chains
%   THIN     - Delete burn-in and thin MCMC-chains
%   JOIN     - Join similar structures of arrays to one structure of arrays
%   BATCH    - Batch MCMC sample chain and evaluate mean/median of batches
%
% Misc:
%   CUSTATS   - Calculate cumulative statistics of data
%   BBPRCTILE - Bayesian bootstrap percentile
%   GRADCHEK  - Checks a user-defined gradient function using finite differences.
% PROBABILITY DISTRIBUTION FUNCTIONS (contents of the dist-folder):
%
% Priors 
%  PRIOR_GAMMA     Gamma prior structure     
%  PRIOR_INVGAM     Inverse-gamma prior structure     
%  PRIOR_LAPLACE    Laplace (double exponential) prior structure
%  PRIOR_LOGLOGUNIF Uniform prior structure for the log(log(parameter))
%  PRIOR_LOGUNIF    Uniform prior structure for the logarithm of the parameter
%  PRIOR_NORMAL     Normal prior structure     
%  PRIOR_SINVCHI2   Scaled inverse-chi-square prior structure
%  PRIOR_T          Student-t prior structure
%  PRIOR_UNIF     Uniform prior structure     
%
% Probability density functions
%
%    BETA_LPDF     - Beta log-probability density function (lpdf).
%    BETA_PDF      - Beta probability density function (pdf).
%    DIR_LPDF      - Log probability density function of uniform Dirichlet
%                    distribution
%    DIR_PDF       - Probability density function of uniform Dirichlet
%                    distribution
%    GAM_CDF       - Cumulative of Gamma probability density function (cdf).
%    GAM_LPDF      - Log of Gamma probability density function (lpdf).
%    GAM_PDF       - Gamma probability density function (pdf).
%    GEO_LPDF      - Geometric log probability density function (lpdf).
%    INVGAM_LPDF   - Inverse-Gamma log probability density function.
%    INVGAM_PDF    - Inverse-Gamma probability density function.
%    LAPLACE_LPDF  - Laplace log-probability density function (lpdf).
%    LAPLACE_PDF   - Laplace probability density function (pdf).
%    LOGN_LPDF     - Log normal log-probability density function (lpdf)
%    LOGT_LPDF     - Log probability density function (lpdf) for log Student's T
%    MNORM_LPDF    - Multivariate-Normal log-probability density function (lpdf).
%    MNORM_PDF     - Multivariate-Normal log-probability density function (lpdf).
%    NORM_LPDF     - Normal log-probability density function (lpdf).
%    NORM_PDF      - Normal probability density function (pdf).
%    POISS_LPDF    - Poisson log-probability density function.
%    POISS_PDF     - Poisson probability density function.
%    SINVCHI2_LPDF - Scaled inverse-chi log-probability density function.
%    SINVCHI2_PDF  - Scaled inverse-chi probability density function.
%    T_LPDF        - Student's T log-probability density function (lpdf)
%    T_PDF         - Student's T probability density function (pdf)
%
% Random number generators
%
%    CATRAND       - Random matrices from categorical distribution.
%    DIRRAND       - Uniform dirichlet random vectors
%    EXPRAND       - Random matrices from exponential distribution.
%    GAMRAND       - Random matrices from gamma distribution.
%    INTRAND       - Random matrices from uniform integer distribution.
%    INVGAMRAND    - Random matrices from inverse gamma distribution
%    INVGAMRAND1   - Random matrices from inverse gamma distribution
%    INVWISHRND    - Random matrices from inverse Wishart distribution.
%    NORMLTRAND    - Random draws from a left-truncated normal
%                    distribution, with mean = mu, variance = sigma2
%    NORMRTRAND    - Random draws from a right-truncated normal
%                    distribution, with mean = mu, variance = sigma2
%    NORMTRAND     - Random draws from a normal truncated to interval
%    NORMTZRAND    - Random draws from a normal distribution truncated by zero
%    SINVCHI2RAND  - Random matrices from scaled inverse-chi distribution
%    TRAND         - Random numbers from Student's t-distribution
%    UNIFRAND      - Generate unifrom random numberm from interval [A,B]
%    WISHRND       - Random matrices from Wishart distribution.
%
% Others
%    KERNELP       - Kernel density estimator for one dimensional distribution.
%    HAMMERSLEY    - Hammersley quasi-random sequence
% MCMC FUNCTIONS (contents of the /mc-folder):
%
%    BBMEAN        - Bayesian bootstrap mean
%    GIBBS         - Gibbs sampling
%    HMC2          - Hybrid Monte Carlo sampling.
%    HMC2_OPT      - Default options for Hybrid Monte Carlo sampling.
%    HMEAN         - Harmonic mean
%    METROP2       - Markov Chain Monte Carlo sampling with Metropolis algorithm.
%    METROP2_OPT   - Default options for Metropolis sampling.
%    RESAMPDET     - Deterministic resampling
%    RESAMPRES     - Residual resampling
%    RESAMPSIM     - Simple random resampling
%    RESAMPSTR     - Stratified resampling
%    SLS           - Markov Chain Monte Carlo sampling using Slice Sampling
%    SLS_OPT       - Default options for Slice Sampling
%    SLS1MM        - 1-dimensional fast minmax Slice Sampling
%    SLS1MM_OPT    - Default options for SLS1MM_OPT
%    SOFTMAX2      - Softmax transfer function
% MISCELLANEOUS FUNCTIONS (contents of the /misc-folder):)
%
%    CVIT         - Create itr and itst indeces for k-fold-cv
%    CVITR        - Create itr and itst indeces for k-fold-cv with ranndom 
%                   permutation
%    MAPCOLOR     - returns a colormap ranging from blue through gray
%                   to red
%    MAPCOLOR2    - Create a blue-gray-red colormap.
%    M2KML        - Converts GP prediction results to a KML file
%    QUAD_MOMENTS - Calculate the 0th, 1st and 2nd moment of a given (unnormalized) 
%                   probability distribution
%    RANDPICK     - Pick element from x randomly
%                   If x is matrix, pick row from x randomly.
%    STR2FUN      - Compatibility wrapper to str2func
%    SET_PIC      - Set the inducing inputs and blocks for two dimensional input data
%    WMEAN        - weighted mean
% OPTIMIZATION FUNCTIONS (contents of the /optim-folder):
%
%    BSEARCH       - Finds the minimum of a combinatorial function using backward search
%    BSEARCH_OPT   - Default options for backward search
%    FSEARCH       - Finds the minimum of a combinatorial function using forward search
%    FSEARCH_OPT   - Default options for forward search
%    SCGES         - Scaled conjugate gradient optimization with early stopping
%    SCGES_OPT     - Default options for scaled conjugate gradient optimization
%    SCGES         - Scaled conjugate gradient optimization with early stopping (new options structure).
%    SCG2          - Scaled conjugate gradient optimization
%    SCG2_OPT      - Default options for scaled conjugate gradient optimization (scg2) (new options structure).
