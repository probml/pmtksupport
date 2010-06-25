#include <stdio.h>
#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <limits.h>
#include <iomanip>
#include <time.h>
#include <sys/stat.h>
#include "mex.h"
#include "various.h"
#include "matrix.h"
extern mwSize SampleSize;
extern mwSize NumberOfGenes;
extern char *DataFile;
extern double **Data;

void CheckPointer(void* pomwSizeer)
{
   if(pomwSizeer==NULL)
   {
      printf("Memory allocation error.\n");
      exit(1);
   }
   return;
}

 
//computes n choose m
mwSize choose(mwSize m, mwSize n)
{
   if(m > n) return 0;
   if(m == 0) return 1;
   if(m == 1) return n;
   if(m == n) return 1;
   if(m == (n-1)) return n;

   return(choose(m, n-1) + choose(m-1, n-1));
}

double logdexp(double x,double beta)
{
   return(-log(beta)-x/beta);
}

double logdpois(double x,double lambda)
{
   mwSize i;
   double s = -lambda+x*log(lambda);

   for(i=2;i<=x;i++)
   {
      s += log((double)i);
   }

   return(s);
}

double getLogBetaPrior(double m,double x)
{
   double tau = 100;
   return((tau*m-1)*log(x)+(tau*(1-m)-1)*log(1-x));
}

double** ReadData()
{
   mwSize i,j;
   double s;
   double** data;

   data = new double*[SampleSize];
   for(i=0;i<SampleSize;i++)
   {
      data[i] = new double[NumberOfGenes];
   }

   FILE* in = fopen(DataFile,"r");
   if(NULL==in)
   {
      printf("Cannot open data file %s.\n",DataFile);
      exit(1);
   }
   for(i=0;i<SampleSize;i++)
   {
      for(j=0;j<NumberOfGenes;j++)
      {
	 fscanf(in,"%lf",&s);
         data[i][j] = s;
      }
   }
   fclose(in);
   return data; 
}

void FreeData(double** Data)
{
   mwSize i;

   for(i=0;i<SampleSize;i++)
   {
      delete[] Data[i]; 
      Data[i] = NULL;
   }
   delete[] Data;
   Data = NULL;
   return;
}

void NormalizeData(double** data)
{
   mwSize i,j;
   double s;

   for(j=0;j<NumberOfGenes;j++)
   {
      s = 0.0;
      for(i=0;i<SampleSize;i++)
      {
	 s += data[i][j];
      }
      s = s/((double)SampleSize);
      for(i=0;i<SampleSize;i++)
      {
	 data[i][j] -= s;
      }
   }
   return;
}

double getSSD(double** data,mwSize vi,mwSize vj)
{
   mwSize k;
   double s = 0.0;

   for(k=0;k<SampleSize;k++)
   {
      s += data[k][vi]*data[k][vj];     
   }
   return s;
}

double mylgamma(double lambda,mwSize p)
{
   mwSize i;
   double s = 0.0;

   for(i=0;i<p;i++)
   {
      s += lgamma(lambda-(((double)i)/2.0));
   }
   return(s);
}


#if defined(_WIN32)
/* Log of gamma function by asymptotic series */
double lgamma(double x){

double Recxs,sum,term,x10;

const double EulerG = 0.577215664902,
 hlfLn2Pi = 0.918938533205;
const double b[] =
 {0.0, 0.0833333333333, -0.00277777777778,
  0.000793650793651, -0.000595238095238,
  0.000841750841751};
mwSize xInt,k,x11;

if ( x > 10.0 ) /* use asymptotic series */
 {
 Recxs = 1.0/(x*x);
 term = x;
 sum = (x-0.5)*log(x)-x+hlfLn2Pi;
 for ( k = 1; k <= 5; k++ )
  {
  term *= Recxs;
  sum += b[k]*term;
  }
 return  sum;
  }

if ( x > 0 ) /* recurrence to x>10 */
 {
 x -= 1.0;
 x11 = (mwSize)(11.0 - x);
 x10 = x + (double) x11;
 xInt = (mwSize)x;
 sum = 0.0;
 for ( k = 1; k <= x11-1; k++ )
        {
        sum -= log(x10 - (double) k );
        }
 return sum+lgamma(x10);
 }
}

double drand48(){
  double x;
  x = ((double)rand())/((double)RAND_MAX);
  return x;   
}
#endif


bool FileExists(char* strFilename) {
  struct stat stFileInfo;
  bool blnReturn;
  mwSize mwSizeStat;

  // Attempt to get the file attributes
  mwSizeStat = stat(strFilename,&stFileInfo);
  if(mwSizeStat == 0) {
    // We were able to get the file attributes
    // so the file obviously exists.
    blnReturn = true;
  } else {
    blnReturn = false;
  }
  
  return(blnReturn);
}
