Author: Jarno Vanhatalo <Jarno.Vanhatalo@tkk.fi>
Last modified: 2010-07-17 08:51:11 EEST
-----------------------------------------------------------------

GPstuff: Gaussian process models for Bayesian analysis V2.0 (r395)

Maintainer: Jarno Vanhatalo <jarno.vanhatalo@tkk.fi>

This software is distributed under the GNU General Public Licence
(version 2 or later); please refer to the file Licence.txt,
included with the software, for details.


Table of contents:

1. INTRODUCTION
2. DOCUMENTATION AND DEMOS
3. INSTALLING THE TOOLBOX
4. USER QUIDE (VERY SHORT)
5. KNOWN PROBLEMS WITH INSTALLING SUITESPARSE AND SOLUTIONS TO THEM


    ------------------------------------------
1. INTRODUCTION

  GPstuff is a collection of Matlab functions to build and analyze
  Bayesian models build over Gaussian processes. The toolbox is tested
  with Matlab 7.9 in 64bit Windows and Linux environments (it should
  work in the 32bit versions as well but they are not tested
  properly).

2. INSTALLING THE TOOLBOX

  Some of the functions in GPstuff are implemented using C in order to
  make the computations faster. In order to use these functions you
  need to compile them first. There are two ways to do that:

  1) Basic installation without compactly supported covariance
     functions

  * Install the GPstuff package by running 'matlab_install' in this
    folder

  * With this option you are able to use all the other functions
    except for gpcf_ppcs*


  2) Installation with compactly supported covariance functions
  
  Compactly supported (CS) covariance functions are functions that
  produce sparse covariance matrices (matrices with zero elements). To
  use these functions (gpcf_ppcs*) you need the sparse GP
  functionalities in GPstuff which are build over SuiteSparse toolbox
  by Tim Davis. To take full advantage of the CS covariance functions
  install GPstuff as follows:

  * First install SuiteSparse from:
    http://www.cise.ufl.edu/research/sparse/SuiteSparse/current/SuiteSparse/

     Note!  Install also Metis 4.0.1 as mentioned in the site under
            header "Other packages required:".
     Note2! There are problems with installing SuiteSparse in some
            architechtures. See the end of this document for known
            problems and solutions.

  * Install the GPstuff package:

    Run 'matlab_install( suitesparse_path )' in the present directory. 
    Here suitesparse_path is a string telling the path to SuiteSparse 
    package (for example, '/matlab/toolbox/SuiteSparse/'). Note! It is
    important that suitesparse_path is in right format. Include also
    the last '/' sign in it. In windows replace / by \.

    The function matlab_install compiles the mex-files and prints on
    the screen, which directories should be added to Matlab paths. 
    

3. CONTENTS
   
   The packge contains the following subdirectories:
   diag  dist  gp  mc  misc  optim

   Each folder contains Contents.m, which summarizes the functions
   in the folder. 

   From the above 'gp' folder contains the main functionalities and
   demonstration programs. Other folders contain additional help
   functions.

4. USER QUIDE (VERY SHORT)

   It easiest to learn to use the package by running the demos. It is
   advisable to open the demo files in text editor and run them line
   by line. The demos are documented so that user can follow what
   happens on each line.

   The basic structure of the program is as follows. The program
   consist of separate blocks, which are:

      Gaussian process model structure (GP):
                      This is a structure that contains all the
                      model information (see gp_init) and information
                      on, which inference scheme is used. 

                      GP structure contains covariance function
                      structures (GPCF_*) and likelihood structures
                      (LIKELIH_*). 

      Covariance function structure (GPCF):
                      This is a structure that contains all of the
                      covariance function information (see
                      e.g. gpcf_sexp). The structure contains the
                      hyperparameter values, pointers to nested
                      functions that are related to the covariance
                      function (e.g. function to evaluate covariance
                      matrix) and hyperprior structure.

      likelihood structure:
                      This is a structure that contains all of the
                      likelihood function information (see
                      e.g. likelih_probit). The structure contains the
                      likelihood parameter values and pointers to
                      nested functions that are related to the
                      likelihood function (e.g. log likelihood and its
                      derivatives).

      Inference utilities:
                      Inference utilities consist of functions that
                      are needed to make the posterior inference and
                      predictions. These include, among others,
		        GP_E - Evaluate conditional log posterior
 		               density
                        GP_G - Evaluate gradient of conditional log
                               posterior 
			EP_PRED - Predictions with Gaussian Process EP
                        GP_MC - Markov chain Monte Carlo sampling



5. KNOWN PROBLEMS WITH INSTALLING SUITESPARSE AND SOLUTIONS TO THEM

  There are some problems with installing SuiteSparse for  64bit Linux,
  Matlab 7.8 (or newer). The problem can be fixed as follows


  - Compiling SuiteSparse seems to finish allright. However, when running
    demos you will end up in segmentation fault. 
  - The problem occurs with CHOLMOD, UMFPACK, and SPQR packages and can
    be prevented with the following steps:
 
  * Open file SuiteSparse/CHOLMOD/MATLAB/cholmod_make.m
  * Go to the line 46
  * Replace the following lines 

    (46) if (~isempty (strfind (computer, '64')))
    (47)     % 64-bit MATLAB
    (48)     d = '-largeArrayDims' ;
    (49) end

    with the following ones

    if (~isempty (strfind (computer, '64')))
        % 64-bit MATLAB
        d = '-largeArrayDims' ;
        
        % The next three lines are added by jarno.vanhatalo@iki.fi (-09).
        % These options are needed for some reason in Matlab 7.8 or newer. 
        if v >= 7.8 
            d = [d ' -DLONG -D''LONGBLAS=UF_long''']; 
        end
    end


  * Open file SuiteSparse/SPQR/MATLAB/spqr_make.m
  * Go to the line 58
  * Replace the following lines

   (58) if (is64)
   (59)    % 64-bit MATLAB
   (60)    d = '-largeArrayDims' ;
   (61) end

   with the following ones 

   if (is64)
       % 64-bit MATLAB
       d = '-largeArrayDims' ;
    
       % The next three lines are added by jarno.vanhatalo@iki.fi (-09).
       % These options are needed for some reason in Matlab 7.8 or newer.
       if v >= 7.8 
           d = [d ' -DLONG -D''LONGBLAS=UF_long''']; 
       end
   end

  * Open file SuiteSparse/UMFPACK/MATLAB/umfpack_make.m
  * Go to the line 21
  * replace the following lines

   (21) if (~isempty (strfind (computer, '64')))
   (22)    d = ' -largeArrayDims' ;
   (23) end

   with the following ones

   v = getversion ;      % Added by jarno.vanhatalo@iki.fi (-09).
   if (~isempty (strfind (computer, '64')))
       d = ' -largeArrayDims' ;
    
        % The next three lines are added by jarno.vanhatalo@iki.fi (-09).
        % These options are needed for some reason in Matlab 7.8 or newer. 
        if v >= 7.8 
            d = [d ' -DLONG -D''LONGBLAS=UF_long''']; 
        end
   end