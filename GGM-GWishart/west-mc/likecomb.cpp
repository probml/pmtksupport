 /* 4-1-03
The purpose of these functions is to compute the likelihood of non-decomposable
graphs by combining the likelihood calculation for cliques and separators, and
the Importance sampling algorithm for prime components that are not complete.
This function will take a list of the prime components, determine which are 
complete, execute the appropriate computation, and then recombine them. 
*/
#include <stdio.h>
#include <stdlib.h>
#include <cstdlib>
#include <math.h>
#include <string.h>
#include <limits.h>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <time.h>
#include <mex.h>
#include "graph.h"
#include "various.h"
#include "ISampling.h"
#include "mex.h"
#include "matrix.h"
using namespace std;

#define PI 3.141592653589793

#if defined(_WIN32)
extern "C" {
void dpotrf(char *, mwSize *, double *, mwSize *, mwSize *);
void dgetrf(mwSize *, mwSize *, double *, mwSize *, mwSize *, mwSize *);
void dgetri(mwSize *, double *, mwSize *, mwSize *, double *, mwSize *, mwSize *);
}
#elif defined(_WIN64)
extern "C" {
void dpotrf(char *, mwSize *, double *, mwSize *, mwSize *);
void dgetrf(mwSize *, mwSize *, double *, mwSize *, mwSize *, mwSize *);
void dgetri(mwSize *, double *, mwSize *, mwSize *, double *, mwSize *, mwSize *);
}
#else
extern "C" {
void dpotrf_(char *, mwSize *, double *, mwSize *, mwSize *);
void dgetrf_(mwSize *, mwSize *, double *, mwSize *, mwSize *, mwSize *);
void dgetri_(mwSize *, double *, mwSize *, mwSize *, double *, mwSize *, mwSize *);
}
#endif


/*this says cliques, but it can be used for separators too*/

double cliquelike(int *vertices, int nvertices, double ** suffdat, double tau, double delta){
  extern double **Data;
  extern mwSize SampleSize;    //was mwSize
  extern mwSize NumberOfGenes; //was mwSize
  extern int debug_flag;
  mwSize mw_nvertices = nvertices; 

  double * Phi_star;
  char UPLO='L';
  mwSize INFO;
  
  double answer;
  mwSize i,j;

   Phi_star=new double[nvertices*nvertices];
  
  answer=0;

  if(nvertices >1){
     /*gamma part of both prior and posterior*/
    for(i=0; i< nvertices; i++){
      answer=answer-lgamma((delta+((double) nvertices-i)-1.0)/2.0);
      answer=answer+lgamma((delta+((double) SampleSize+nvertices-i)-1.0)/2.0);
      }
    /*Determinant of Phi */
    answer=answer+log(tau)*((double) nvertices)*(delta+nvertices-1.0)/2.0;

      
    mwSize index=0;
    for(i=0; i< nvertices; i++){
      for(j=0; j< nvertices; j++){
        Phi_star[index]=suffdat[vertices[i]][vertices[j]];
        if(i==j) Phi_star[index]=Phi_star[index]+tau;
        index++;
      }
    }

    #if defined(_WIN32)    
    dpotrf(&UPLO, &mw_nvertices, Phi_star, &mw_nvertices, &INFO);
    #elif defined(_WIN64)    
    dpotrf(&UPLO, &mw_nvertices, Phi_star, &mw_nvertices, &INFO);
    #else
    dpotrf_(&UPLO, &mw_nvertices, Phi_star, &mw_nvertices, &INFO);
    #endif
                  
    /*determinant of Phi_Star*/     
      for(i=0; i< nvertices; i++){
        index=i*nvertices+i;
        answer=answer-log(Phi_star[index])*(delta+((double) nvertices+SampleSize)-1.0);
       }      
  }

  else{
    answer += lgamma(((double)(delta+SampleSize))/2.0);
    answer -= lgamma(delta/2.0);
    answer += delta*log(tau)/2.0;
    answer -= ((double)(delta+SampleSize))*log(tau+suffdat[vertices[0]][vertices[0]])/2.0;
  }


 
  delete[] Phi_star;


  return answer;
   
}



int complete(int ** incidence, int * vertices, int nvertices){
int i,j, answer;
answer=1;
 if(nvertices >1){
   i=0;	
   while(answer==1 && i<nvertices){
     j=i+1;
     while(answer==1 && j<nvertices){
       if(incidence[vertices[i]][vertices[j]]!=1) answer=0;
       j++;}
     i++;
   }}

return answer;
}

double likecomb(Graph *graph, double ** suffdat, double tau, double delta){
  extern double **Data;
  extern int debug_flag;
  extern mwSize nsamples, ntrials;
  extern mwSize SampleSize;
  extern mwSize NumberOfGenes;
  mwSize i,j;
  double answer, posNC, priNC;

  answer=-((double)SampleSize)*((double)NumberOfGenes)*log(PI)/2.0;

  for(i=0; i<graph->nCliques; i++){
   
   if(graph->CliquesDimens[i]< SampleSize){
     if( complete(graph->Edge, graph->Cliques[i], graph->CliquesDimens[i])){  
       posNC=cliquelike(graph->Cliques[i], graph->CliquesDimens[i], suffdat, tau, 
	   delta);
       answer=answer+posNC;
     }
     else{  
       mwSize t,k, posconverged = 0,priconverged =0;
       double precision = 1e-2;
       double *pvals = new double[ntrials];
       double pmean=0;
       double pstd=0;
       double posNC;
       double priNC; 
       char buff[100];
       string fname;

       //Run up to ntrials sampling runs of nsmaples each. 
       //Stop when the standard error of the mean of the estimates 
       //is less than precision*abs(pmean) giving approximately 
       //1/precision significant figures
       pmean = 0;
       for(t=0;t<ntrials;t++){ 
          pvals[t]=Isampling(suffdat,graph->CliquesDimens[i], graph->Cliques[i], 
                  delta+SampleSize, tau, graph->Edge, nsamples);
          pmean = (t*pmean + pvals[t])/(double(t+1));
          pstd  = 0;
          for(k=0;k<=t;k++){
            pstd = pstd + (pvals[k]-pmean)*(pvals[k]-pmean);
          } 
          pstd = sqrt(pstd/double(t+1))/sqrt(double(t+1)); 
          if((t>0) & (pstd<=precision*fabs(pmean))){
            posconverged = 1;
            break;
          }
       }
       posNC = pmean;
       if(!posconverged & debug_flag){
         mexPrintf("Warning: Numerator did not converge in %d trials: %.4e %.4e\n", ntrials,posNC, pstd);
       }
  
       pmean = 0;
       for(t=0;t<ntrials;t++){ 
          pvals[t]=Isampling(NULL, graph->CliquesDimens[i], graph->Cliques[i], 
                      delta, tau, graph->Edge, nsamples);
          pmean = (t*pmean + pvals[t])/(double(t+1));
          pstd  = 0;
          for(k=0;k<=t;k++){
            pstd = pstd + (pvals[k]-pmean)*(pvals[k]-pmean);
          } 
          pstd = sqrt(pstd/double(t+1))/sqrt(double(t+1)); 
          if(t>0 & pstd<=precision*fabs(pmean)){
            priconverged = 1;
            break;
          }
       }
       priNC = pmean;
       if(!priconverged & debug_flag){
        mexPrintf("Warning: Denominator did not converge in %d trials: %.4e %.4e\n", ntrials,priNC, pstd);
       }
         
  
        if((!priconverged | !posconverged) & debug_flag){
          mexPrintf("Warning: did not converge in %d trials: pri:%.4e pos:%.4e\n", ntrials,priNC, posNC);
        }

       answer=answer+posNC-priNC;

       delete[] pvals;
     }
   }
else{
 answer= -100000.00;
 break;}
 }

for(i=0; i< graph->nSeparators; i++){
posNC=cliquelike(graph->Separators[i], graph->SeparatorsDimens[i], suffdat, tau, delta);
 answer=answer-posNC;
 }

return answer;

}



