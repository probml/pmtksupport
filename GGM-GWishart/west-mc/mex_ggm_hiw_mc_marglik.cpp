//using namespace std;
 
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <limits.h>
#include <iostream>
#include <iomanip>
#include <time.h>
#include "graph.h"
#include "likecomb.h"
#include "various.h"
#include "graphfns.h"
#include "ISampling.h"
#include "mex.h"
#include "matrix.h"

int debug_flag  = 0;     //Set to 1 to enable MC convergence messages
mwSize ntrials;          //Max number of sampling trials to reach convergence
mwSize nsamples;         //Number of samples drawn per sampling trial

mwSize NumberOfGenes;  //columns in dataset (was INT)
mwSize SampleSize;     //rows in dataset (was INT)
double **Data=NULL;
char *DataFile;

void mexFunction(int nlhs, mxArray *plhs[ ],int nrhs, const mxArray *prhs[ ]) {

  double *marglik;
  mwSize i, j, k, total_edges, current_edges, add;	
  mwSize seed=0;	 
  double tau;
  mwSize * edge;
  double delta;
  double beta =0;
  char *GraphFile;
  mwSize n_iterations = 0;
  double *tmp;
  int *tmpi;
  double *cov;
  int *edgelist;
  mwSize numedge;

  /* Check for correct number of input/output arguments*/
  if (nrhs != 8) {
	mexErrMsgTxt("Eight input arguments required.");
  } 
  if(nlhs > 1){
	mexErrMsgTxt("Too many output arguments.");
  }
  plhs[0] = mxCreateDoubleMatrix(1,1,mxREAL);
  marglik = mxGetPr(plhs[0]);

  /* Check data type of input argument  1*/
  if (!(mxIsDouble(prhs[0]))){
    mexErrMsgTxt("Input argument scattermatrix must be of type double.");
  }
  cov = mxGetPr(prhs[0]); 
  NumberOfGenes = (mwSize)mxGetM(prhs[0]);

  /* Check data type of input argument  2*/
  if (!(mxIsInt32(prhs[1]))){
    mexErrMsgTxt("Input argument graph must be of type int32.");
  }
  edgelist = (int*)mxGetData(prhs[1]);
  numedge  = (mwSize)mxGetM(prhs[1]); 

  /* Check data type of input argument  3*/
  if (!(mxIsDouble(prhs[2]))){
    mexErrMsgTxt("Input argument samplesize must be of type double.");
  }
  tmp = mxGetPr(prhs[2]);
  SampleSize = (mwSize)*tmp;

  /* Check data type of input argument  4*/
  if (!(mxIsDouble(prhs[3]))){
    mexErrMsgTxt("Input argument tau must be of type double.");
  }
  tmp = mxGetPr(prhs[3]);
  tau = *tmp;

  /* Check data type of input argument  5*/
  if (!(mxIsDouble(prhs[4]))){
    mexErrMsgTxt("Input argument delta must be of type double.");
  }
  tmp = mxGetPr(prhs[4]);
  delta = *tmp;

  /* Check data type of input argument  6*/
  if (!(mxIsInt32(prhs[5]))){
    mexErrMsgTxt("Input argument ntrials must be of type int32.");
  }
  tmpi = (int*)mxGetData(prhs[5]); 
  ntrials = *tmpi;
  
  /* Check data type of input argument  7*/
  if (!(mxIsInt32(prhs[6]))){
    mexErrMsgTxt("Input argument nsamples must be of type int32.");
  }
  tmpi = (int*)mxGetData(prhs[6]); 
  nsamples = *tmpi; 
  
  /* Check data type of input argument  8*/
  if (!(mxIsUint32(prhs[7]))){
    mexErrMsgTxt("Input argument randseed must be of type uint32.");
  }
  tmpi = (int*)mxGetData(prhs[5]);
  #if defined(_WIN32)
  srand(*tmpi);
  #else
  srand48(*tmpi);
  #endif

  double current_loglike, proposed_loglike,hastings_ratio, cutoff;
  double **suffdat = new double*[NumberOfGenes];
  

  /*Copy sufficient stats from flat matlab array to 2D C array*/
  for(i=0; i< NumberOfGenes; i++){
    suffdat[i]=new double[NumberOfGenes];
    for(j=0; j< NumberOfGenes; j++){
      suffdat[i][j]=cov[i + j*NumberOfGenes];
    }
  }

  /*Initialize the graph structure*/
  edge=new mwSize[2];
  LPGraph graph = new Graph;
  total_edges=(NumberOfGenes*(NumberOfGenes-1))/2;

  /*Create the graph structure*/ 
  graph->EdgeListToMss(edgelist,(int)numedge);
  graph->InitGraphFromMss();
  graph->GetMPSubgraphs();
  graph->InitConnectedComponents();
  current_edges=0;
  for(i=0; i< NumberOfGenes; i++){
    for(j=(i+1); j< NumberOfGenes; j++){
      if(graph->Edge[i][j]==1) current_edges ++;
    }
  }

  //Compute the marginal likelihood estimate 
  *marglik =likecomb(graph, suffdat, tau, delta);

  //Free 2D suffdat array
  for(i=0; i< NumberOfGenes; i++){
    delete [] suffdat[i];
  }
  delete [] suffdat;
  delete [] edge;

  return;
}






