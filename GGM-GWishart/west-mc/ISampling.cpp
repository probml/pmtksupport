#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <limits.h>
#include <iostream>
#include "mex.h"
#include "random.h"
#include "various.h"
#include <float.h>
#ifndef _GRAPH_H
#include "graph.h"
#include "mex.h"
#include "matrix.h"
#endif
#define PI 3.141592653589793
using namespace std;

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




double Isampling(double** suffdata, mwSize size, int* vert, double delta, double tau, int** A, mwSize N)
{
   mwSize ii;   
   extern mwSize NumberOfGenes;
   extern int debug_flag;
   mwSize nvar=NumberOfGenes;
   double *T = new double[size*size];
   double *tt = new double[size*size];
   double *data = new double[size*size];
   mwSize *AA = new mwSize[size*size];
   mwSize *nu = new mwSize[size];
   mwSize *kk = new mwSize[size];
   mwSize *b  = new mwSize[size];

   mwSize i,j,k,s,ss;
   char *UPLO=new char[1];
   UPLO[0]='L'; 
   mwSize *IPIV;
   mwSize INFO=0;  
   mwSize lwork=3*size;
   double max_is, min_is, trial_run[1000], offset,eesq;
   double *work = new double[lwork];


   IPIV= new mwSize[size];
   for(i=0;i<size;i++){
     IPIV[i] = 0;
   }

   for(i=0;i<size;i++)
   {  
      j=vert[i];
      for(k=0;k<size;k++)
      {
         s=vert[k];
	 if(suffdata!=NULL) data[i*size+k]=suffdata[j][s];
	 else data[i*size+k]=0.0;
	 if(i==k) data[i*size+k] = data[i*size+k]+tau;
      }
   }
    
   for(i=0;i<size;i++)
   {
      j=vert[i];
      for(k=0;k<size;k++)
      {
         s=vert[k];
         AA[i*size+k]=A[j][s];
	 if(i==k) AA[i*size+k]=1;
      }
   } 
    
   for(i=0;i<size;i++)
   {
      nu[i]=0;
      for(j=(i+1);j<size;j++)
      {
         nu[i]+=AA[i*size+j];
      }
   }

   //  for(i=0;i<size;i++){cout<<nu[i]<<" ";}
    
   for(i=0;i<size;i++)
     {	kk[i]=0;
      for(j=0;j<i;j++)
      {
         kk[i]+=AA[j*size+i];
	 }
   }

   //for(i=0;i<size;i++){cout<<kk[i]<<" ";}
   
   for(i=0; i<size; i++) 
   {
      for(j=0; j<size; j++) 
      {
         T[i*size+j] = data[j*size+i];
      }
   }

   //DPOTRF( UPLO, N, A, LDA, INFO ) Cholesky decomposition from LAPACK!!!
   //DGETRF( M, N, A, LDA, IPIV, INFO ) LU Factorization
   //DGETRI( N, A, LDA, IPIV, WORK, LWORK, INFO ) Inverse...       
   
    #if defined(_WIN32)
    dgetrf(&size,&size,T,&size,&IPIV[0],&INFO);
    dgetri(&size,T,&size,&IPIV[0],work,&lwork,&INFO);
    dpotrf(UPLO,&size,T,&size,&INFO);
    #elif defined(_WIN64)
    dgetrf(&size,&size,T,&size,&IPIV[0],&INFO);
    dgetri(&size,T,&size,&IPIV[0],work,&lwork,&INFO);
    dpotrf(UPLO,&size,T,&size,&INFO);    
    #else
    dgetrf_(&size,&size,T,&size,&IPIV[0],&INFO);
    dgetri_(&size,T,&size,&IPIV[0],work,&lwork,&INFO);
    dpotrf_(UPLO,&size,T,&size,&INFO);
    #endif
        
            
   for(i=1;i<size;i++)
   {
      for(j=0;j<i;j++)
      {
         T[i*size+j]=0;
      }
   }

    
   for(i=0;i<size;i++)
   {
      for(j=0;j<size;j++)
      {
         tt[i*size+j]=T[i*size+j]/T[j*size+j];
      }
   }

   for(i=1;i<size;i++)
   {
      for(j=0;j<i;j++)
      {
         tt[i*size+j]=0;
      }
   }

   
   for(i=0;i<size;i++)
   {
      b[i]=nu[i]+kk[i]+1;
   }
   
  double ee=0,lixo,aux,aux1,aux2,is; 
  double *phi = new double[size*size];
  
  eesq=0;
  double msqee;
  mwSize enough=0, oldN=0;


  for(i=0;i<N;i++)
  {
     for(j=0;j<size;j++)
     {
        for(k=0;k<size;k++)
        {
           phi[j*size+k]=0;
	}
     }
  
     
     //STEP 3

     for(j=0;j<size;j++)
     { 
        phi[j*size+j]=sqrt(rchisq(delta+((double) nu[j]))); 
     }
    
     for(j=0;j<(size-1);j++)
     {
        for(k=(j+1);k<size;k++)
        {
           if(AA[j*size+k]!=0)
           { 
              phi[j*size+k]=rnor(0.0,1.0); 
	   }
	}
     }

     //STEP 4

     for(j=1;j<size;j++)
     {
        if(AA[0*size +j]==0)
        {
           for(k=0;k<j;k++)
           {
	     phi[0*size+j]+=phi[0*size + k]*tt[k*size+j];   //OK, starts at zero.
	   }  
           phi[0*size+j]= -phi[0*size+j];     
	                                     
	}
     }

     

     for(j=1;j<(size-1);j++)
     {
        for(k=(j+1);k<size;k++)
        {
           

           if(AA[j*size+k]==0)
           {
         
	     lixo=0;
              for(s=0;s<j;s++)
	      {  aux=0;
	      for(ii=s;ii<j;ii++)  
                 {
                    aux+=phi[s*size+ii]*tt[ii*size+j];
	         }
                 aux1=0;
                 //lixo=0;   //shoud be set before 1st loop?
                 for(ii=s;ii<k;ii++)
                 {
                    aux1+=phi[s*size+ii]*tt[ii*size+k];
		   
	         }
	  
	         lixo+=((phi[s*size+j]+aux)/phi[j*size+j])*(phi[s*size+k]+aux1);
		 

	      }
	      aux2=0;
              for(ii=j;ii<k;ii++)  
              {
                 aux2+=phi[j*size + ii]*tt[ii*size + k];
		
	      }
	      
              phi[j*size+k]=-aux2-lixo;
	      

	   }
        }
     }

    
     
     //STEP 5


    
     is=0;

     for(j=0;j<(size-1);j++)
     {
       for(k=(j+1);k<size;k++)          
	 { 
           if(AA[j*size+k]!=1)
           {
	     if(phi[j*size+j]< sqrt(DBL_MAX-is))
	     is+=phi[j*size+k]*phi[j*size+k];  
	     else{
	       is=DBL_MAX;
	       break;}
	   }
        }
       if(is==DBL_MAX) break;
     }            
     if(i < 1000){
       trial_run[i]=is;
       if(i==0){
	 max_is=is;
	 min_is=is;}
       else if(is > max_is) max_is=is;
       else if (is < min_is) min_is=is;	 
       if(i==999){
	 offset=min_is;
	 for(j=0 ;j <=999 ;j++)
	   if(trial_run[j]-offset < 1400.0)
	     ee+=exp((-0.5)*(trial_run[j]-offset));
       } } 
     else{  
       if(is-offset < 1400.0)
	 ee+=exp((-0.5)*(is-offset));  //will force ee to be at least exp(-.5)
     }
  }
  ee=ee/(N);

  
  double CC;
  CC=0;
  for(i=0;i<size;i++)
  {
    
    CC+=(lgamma((delta+nu[i])/2.0))+(log(T[i*size+i])*(delta+b[i]-1));
        
  }
   
  
  if(ee==0 & debug_flag){ 
    mexPrintf("warning, 0  wt. min_is=%.10f \n", min_is);
  }

  ee=log(ee)-offset/2.0;


  delete [] T;
  delete [] tt;
  delete [] data;
  delete [] AA;
  delete [] nu;
  delete [] kk;
  delete [] b;
  delete [] work;
  delete [] phi;
  delete [] IPIV;

  return CC+ee;

}
   								   











