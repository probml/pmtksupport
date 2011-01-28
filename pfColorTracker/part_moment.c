/*

  Compute empirical mean & covariance of (weighted) data

  
  usage: [Smean , Pcov] = part_moment(S , w);
  -----

  Inputs
  ------

  S             Data (d x N)
  w             weights (1 x N)


  Ouputs
  -------

  Smean         Mean vector (d x 1)
  Pcov          Covariance matrix (d x d)


   To compile
  -----------
 
	  
mex  -output part_moment.dll part_moment.c

mex -f mexopts_intel10.bat -output part_moment.dll part_moment.c


  Author : Sébastien PARIS  © (sebastien.paris@lsis.org) (2004)
  -------

  
*/

#include "mex.h"
#include "math.h"


#ifndef max
    #define max(a,b) (a >= b ? a : b)
    #define min(a,b) (a <= b ? a : b)
#endif

/*-------------------------------------------------------------------------------------------------*/


void part_moment(double * , double * , double * , double * , double * , int , int  , int  );


/*-------------------------------------------------------------------------------------------------*/


void mexFunction(int nlhs, mxArray *plhs[],
				 int nrhs,  const mxArray *prhs[])
{
	
	
	double *S , *w;
	
	double *Smean , *Pcov;
	
	double *Stemp;
	
	
	const int *dimsS, *dimsw;
	
	int *dimsSmean , *dimsPcov;
	
	
	int numdimsS , numdimsw;
	
	int numdimsSmean , numdimsPcov;
	
	int i;
	
	int  d , N , slice=1;
	
	
	
	
	if((nrhs < 1) ||  (nrhs > 3)  )    
		
	{
		mexErrMsgTxt("Usage: [Smean , Pcov] = part_moment(S , w);");
		
	}
	/* --- Input 1 ---*/
	
    S        = mxGetPr(prhs[0]);
    
    numdimsS = mxGetNumberOfDimensions(prhs[0]);
    
	dimsS    = mxGetDimensions(prhs[0]);
	
	
	
	d        = dimsS[0];
	
	N        = dimsS[1];
	
	
	for (i = 2 ; i < numdimsS ; i++)
		
	{
		
		slice *= dimsS[i];
		
		
	}
	
	
	/* --- Input 2 ---*/
	
    w        = mxGetPr(prhs[1]);
    
    numdimsw = mxGetNumberOfDimensions(prhs[1]);
    
	dimsw    = mxGetDimensions(prhs[1]);
	
	if (numdimsw != 2)
		
	{	
		
		mexErrMsgTxt("Second input be (1 x N) or (N x 1");
		
	}
	
	if (max(dimsw[0] , dimsw[1]) != N)
		
	{	
		
		mexErrMsgTxt("Second input be (1 x N) or (N x 1");
		
	}
		
	
	/* --- Output 1 ---*/
	
	
    dimsSmean  = (int *)mxMalloc(numdimsS*sizeof(int));
	
	dimsSmean[0] = d;
	
	dimsSmean[1] = 1;
	
	for (i = 2 ; i < numdimsS ; i++)
		
	{
		
		dimsSmean[i] = dimsS[i];
		
	}
	
	numdimsSmean   = numdimsS;
	
	
	plhs[0]        = mxCreateNumericArray(numdimsSmean, dimsSmean, mxDOUBLE_CLASS, mxREAL);
	
	Smean          = mxGetPr(plhs[0]);
	
	
    dimsPcov       = (int *)mxMalloc(numdimsS*sizeof(int));
	
	dimsPcov[0]    = d;
	
	dimsPcov[1]    = d;
	
	for (i = 2 ; i < numdimsS ; i++)
		
	{
		
		dimsPcov[i] = dimsS[i];
		
	}
	
	numdimsPcov    = numdimsS;
	
	
	plhs[1]        = mxCreateNumericArray(numdimsPcov, dimsPcov, mxDOUBLE_CLASS, mxREAL);
	
	Pcov           = mxGetPr(plhs[1]);
	
	
	
	
	Stemp          =  (double *)mxMalloc(slice*d*N*sizeof(double));
	
	
	
	
	/*---------- Main Call --------- */
	
    
	part_moment(S , w , Smean ,  Pcov , Stemp , d , N , slice );
	
	/*------------------------------ */
	
	
	mxFree(Stemp);
	
	mxFree(dimsSmean);
	
	mxFree(dimsPcov);
	
}

/*-------------------------------------------------------------------------------------------------*/



void part_moment(double *S , double *w , double *Smean , double *Pcov , double *Stemp  , int d , int N , int slice )

{
	
	
	
	
	int i , j , t ,  v;
	
	int jd , vd , vdd , td , vdN , id , jdvdd , jvdd;
	
	double temp;
	
	
	for (v = 0 ; v < slice ; v++)
		
	{
		
		
		vd  = v*d;
		
		vdd = vd*d;
		
		vdN = vd*N;
		
		for (j = 0 ; j < d ; j ++)
			
		{
			
			jd   = j + vdN;
			
			temp = 0.0;
			
			
			for (i = 0 ; i < N ; i++)
				
			{
				
				temp += S[jd + i*d]*w[i]; 
				
			}
			
			
			Smean[j + vd] = temp;
			
			
			for (i = 0 ; i < N ; i++)
				
			{
				id              = jd + i*d;
				
				Stemp[id]       = S[id] - temp ; 			
				
			}
			
		}
		
		for (j = 0 ; j < d ; j ++)
			
		{
			
			jd   = j + vdN;
			
			for (t = 0 ; t <= j ; t++)
				
			{
				temp = 0.0;
				
				td   = t + vdN;
				
				
				for (i = 0 ; i < N ; i++)
					
				{
					id    = i*d;
					
					temp += w[i]*Stemp[jd + id]*Stemp[td + id];
					
				}
				
				Pcov[j + t*d + vdd] = temp;
				
			}		
			
		}
		
		for (j = 0 ; j < d - 1 ; j++)
			
		{
			
			jdvdd  = j*d + vdd;
			
			jvdd   = j   + vdd; 
			
			for (t = j + 1 ; t < d ; t++)
				
			{
				
				Pcov[t*d + jvdd] = Pcov[t + jdvdd];			
				
			}
			
		}
		
	}
	
}
