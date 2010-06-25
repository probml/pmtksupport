Version 1.0 - June 24, 2010

Description: This code package implements several marginal likelihood
approximation methods for non-decomposable Gaussian Graphical Models (GGMs) 
under a G-Wishart prior as described in [1]. The methods include BIC score, 
Laplace approximation (related to [3]), diagonal Hessian Laplace approximation, 
and the Monte Carlo estimator of Atay-Kayis & Massam [2]. 

Package Contents: The main files in included in this package are:

* GWExample.m: runs each method on the Iris Virginica example described in 
Atay-Kayis & Massam [2, p. 333, Figure 3].

* GWScore.m: common wrapper function that computes the score (marginal
likelihood approximation) using any of the four methods.

* GWFit.m: code for finding the posterior mode of the precision matrix under
a G-Wishart prior. Used by BIC and Laplace methods.

Detailed comments about the algorithms and sources for each of the methods are 
contained in the source code.

System Requirements: 

* This package requires Matlab. It has been tested using both 32bit and 64bit
versions of Matlab 7.8+ on Linux and Windows platforms. 

* This package depends on Mark Schmidt's freely available minFunc
optimization package. minFunc must be installed and on the Matlab path for
this package to work. See [4] for details on obtaining minFunc.

* The Monte Carlo estimation code is implemented as a C++ mex file. The mex
code has been verified to compile correctly under the gcc and MS Visual
Studio 2008 compilers on Linux and Windows and mex binaries are included. Use
of other platforms and compilers may require changes to the make file. The
c++ mex code is based on code released by Mike West's group at Duke [5]. 

Contact:

* For questions, comments, bug reports, etc... please contact Ben Marlin at
bmarlin[at]cs[dot]ubc[dot]ca.
 
References:

[1] Baback Moghaddam, Benjamin M. Marlin, Emtiyaz Khan, Kevin Murphy.
Accelerating Bayesian Structural Inference for Non-Decomposable Gaussian
Graphical Models. Proceedings of the Neural Information Processing Systems
Conference (NIPS), 2009.

[2] Atay-Kayis & Massam. A Monte Carlo Method for Computing the Marginal
likelihood in Non-decomposable Gaussian Graphical Models. Biometrika. 2005.

[3] Lenkoski & Dobra, Bayesian structural learning & estimation in GGMs,
University of Washington, Dept. of Statistics, TR-545, 2009.

[4] http://www.cs.ubc.ca/~schmidtm/Software/minFunc/minFunc.html 

[5] http://www.stat.duke.edu/research/software/west/ggm.html