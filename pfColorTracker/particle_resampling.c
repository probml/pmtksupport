/* particle_resampling.c 

  usage: index    = particle_resampling(pdf)
  -----
  
  Inputs
  ------
 
  pdf        Empirical probability density function  (vector (1 x N) or (1 x N))
	
  Ouputs
  -------
	  
  index      Index of particle redisdribution (1 x N)
		
		  
  Example 
  -------
  

  N         = 500;
  pdf       = rand(1 , N);
  pdf       = pdf/sum(pdf);
  index     = particle_resampling(pdf);
				
  To compile
  -----------

  mex -DranSHR3 -f mexopts_intel10.bat -output particle_resampling particle_resampling.c


  Author : Sébastien PARIS  © (sebastien.paris@lsis.org) (2004)
  -------

				  
*/




#include <math.h>
#include "mex.h"
#include "time.h"



/*---------------- Basic generators definition ------------------- */

#define mix(a , b , c) \
{ \
	a -= b; a -= c; a ^= (c>>13); \
	b -= c; b -= a; b ^= (a<<8); \
	c -= a; c -= b; c ^= (b>>13); \
	a -= b; a -= c; a ^= (c>>12);  \
	b -= c; b -= a; b ^= (a<<16); \
	c -= a; c -= b; c ^= (b>>5); \
	a -= b; a -= c; a ^= (c>>3);  \
	b -= c; b -= a; b ^= (a<<10); \
	c -= a; c -= b; c ^= (b>>15); \
}

#define znew   (z = 36969*(z&65535) + (z>>16) )
#define wnew   (w = 18000*(w&65535) + (w>>16) )
#define MWC    ((znew<<16) + wnew )
#define SHR3   ( jsr ^= (jsr<<17), jsr ^= (jsr>>13), jsr ^= (jsr<<5) )
#define CONG   (jcong = 69069*jcong + 1234567)
#define KISS   ((MWC^CONG) + SHR3)




#ifdef ranKISS
 #define randint KISS
 #define rand() (randint*2.328306e-10)
#endif 



#ifdef ranSHR3
 #define randint SHR3
 #define rand() (0.5 + (signed)randint*2.328306e-10)
#endif 


#ifndef max
    #define max(a,b) (a >= b ? a : b)
    #define min(a,b) (a <= b ? a : b)
#endif
/*--------------------------------------------------------------- */

#ifdef __x86_64__
    typedef int UL;
#else
    typedef unsigned long UL;
#endif

/*--------------------------------------------------------------- */


static UL jsrseed = 31340134 , jsr;

#ifdef ranKISS

 static UL z=362436069, w=521288629, jcong=380116160;

#endif


 /*--------------------------------------------------------------- */
 
 
 
 
 void randini(void);  
 
 
 
 void particle_resampling(double * , double * , double * , double * , int );
 
 
 
 /*--------------------------------------------------------------- */
 
 
 
 void mexFunction( int nlhs, mxArray *plhs[],
				 int nrhs, const mxArray *prhs[] )
 {
	 
	 
	 double *pdf , *cdf , *uu;
	 
	 double *indice;
	 
	 int  N , M , n;
	 
	 
	 
	 /* Check input */
	 
	 if(nrhs != 1)
		 
	 {
		 
		 mexErrMsgTxt("indice = multinomiale_resampling(p)");
		 
	 }
	 
	 /*  On récupère la taille du vecteurs u_ord */
	 
	 M       = mxGetM(prhs[0]);
	 
	 N       = mxGetN(prhs[0]);
	 
	 if( ((N != 1) & (M > 1) ) |  ((M != 1) & (N > 1) ) )
	 {
		 
		 mexErrMsgTxt("p must be (N x 1) or (1 x N).");
		 
	 }
	 
	 
	 n       = max(M , N);
	 
	 
	 /* Input 1 */
	 
	 
	 pdf     = mxGetPr(prhs[0]);
	 
	 
	 /* Output 1 */
	 
	 
	 
	 plhs[0] = mxCreateDoubleMatrix(M , N ,  mxREAL);
	 
	 indice   = mxGetPr(plhs[0]);
	 
	 
	 /* vecteur temporaire */
	 
	 
	 
	 cdf       = (double *)mxMalloc(n*sizeof(double)); 
	 
	 uu        = (double *)mxMalloc(n*sizeof(double));
	 
	 
	 /* Rand ~U[0,1] Seed initialization */
	 
	 
	 randini();	
	 
	 
	 /* Resampling */
	 
	 
	 particle_resampling(pdf , cdf , uu , indice , n);
	 
	 
	 /* Free ressources */
	 
	 
	 mxFree(uu);
	 
	 mxFree(cdf);
	 
	 
 }
 
 
 
 /* ----------------------------------------------------------------------- */
 
 
 
 void randini(void) 
	 
 {
	 
	 
	 
	 /* SHR3 Seed initialization */
	 
	 jsrseed  = (UL) time( NULL );
	 
	 jsr     ^= jsrseed;
	 
	 
	 
	 /* KISS Seed initialization */
	 
#ifdef ranKISS
	 
	 
	 z        = (UL) time( NULL );
	 
	 w        = (UL) time( NULL ); 
	 
	 jcong    = (UL) time( NULL );
	 
	 mix(z , w , jcong);
	 
#endif 
	 
	 
 }
 
 
 /* ----------------------------------------------------------------------- */
 
 
 void particle_resampling(double *pdf  , double *cdf ,  double *uu , double *indice , int N)	 
 {


	 double sumcdf , sumuu , one = 1.0 , tmp;

	 int    i  , j;


	 /* Compute Scaled empirical CDF from PDF & Sorted Uniform samples*/


	 cdf[0]    = pdf[0];

	 uu[0]     = -log(rand());


	 sumcdf    = cdf[0];

	 sumuu     = uu[0];


	 for(i = 1 ; i < N ; i++)

	 {
		 tmp       = pdf[i] + cdf[i - 1];

		 cdf[i]    = tmp;

		 tmp       = uu[i - 1] - log(rand());

		 uu[i]     = tmp;


		 sumcdf   += cdf[i];

		 sumuu    += uu[i];		 

	 }


	 sumcdf = one/sumcdf;

	 sumuu  = one/sumuu;


	 /* Normalize */


	 for (i = 0 ; i < N ; i++)

	 {
		 cdf[i]            *= sumcdf;

		 uu[i]             *= sumuu;

	 }


	 /*  On recopie les N premiers éléments de u_ord dans le vecteur uu. Le N + 1 élément est égal à 2. */

	 i = 0;

	 j = 0;

	 while (i < N)

	 {
		 if ( (uu[i] < cdf[j]) | (j == N - 1))

		 {
			 indice[i] = j + 1;

			 i++;
		 }

		 else

		 {
			 if (j < N)
			 {
				 j++;
			 }
		 }
	 }	 
 }
