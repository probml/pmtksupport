/*

  
	mex -f mexopts_intel10.bat -output ellipse.dll ellipse.c 
	
	mex -f mexopts_intel10.bat ellipse.c
	  
    Author : Sébastien PARIS
		


 m        = [130 135; 260 270];
 e        = [20 15; 20 10; 0  3];
 [x , y]  = ellipse(m , e);
    
 plot(x , y);
			
			  
*/

#include <math.h>
#include "mex.h"

#define PI 3.14159265358979323846


/*--------------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------------*/

void ellipse(double * , double * , int , double * , double * , int , double * , double *);

/*--------------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------------*/


void mexFunction( int nlhs, mxArray *plhs[] , int nrhs, const mxArray *prhs[] )

{
	
	
	double *m , *e;
	
	
	int  N;
	
	double *x , *y;
	
	double *cos_angle , *sin_angle;
	
	int  *dimsm , *dimse;
	
	int  *dimsx;
	
	int i , d , slice = 1 , numdimsm , numdimse , numdimsx;
	
	
	/*--------------------------------------------------------------------------------*/
	/*--------------------------------------------------------------------------------*/
	/* -------------------------- Parse INPUT  -------------------------------------- */
	/*--------------------------------------------------------------------------------*/
	/*--------------------------------------------------------------------------------*/
	
	if(nrhs < 2)
		
	{     
		mexErrMsgTxt("At least 2 inputs argument are required for ellipse");	
	}
	
	
	/* ----- Input 1 ----- */
	
	m           = mxGetPr(prhs[0]);
	
	numdimsm    = mxGetNumberOfDimensions(prhs[0]);
	
	dimsm       = mxGetDimensions(prhs[0]);
	
	if ( (dimsm[0] != 2) )
		
	{
		mexErrMsgTxt("m must be at least(d x s1 x .... x sp), d >= 2");	
	}
	
	d              = dimsm[0];
	
	for (i = 1 ; i < numdimsm ; i++)
		
	{
		
		slice *= dimsm[i];
		
	}
	
	
    /* ----- Input 2 ----- */
	
	
	e           = mxGetPr(prhs[1]);
	
	numdimse    = mxGetNumberOfDimensions(prhs[1]);
	
	dimse       = mxGetDimensions(prhs[1]);
	
	if ( (dimse[0] != 3)  )
		
	{
		mexErrMsgTxt("e must be (d x s1 x .... x sp)");	
	}
	
	
	/* ----- Input 5 ----- */
	
	if(nrhs < 3)
		
	{
		
		N = 50;
		
	}
	
	else
		
	{
		
		N  = (int) mxGetScalar(prhs[2]);
		
	}
	
	/*--------------------------------------------------------------------------------*/
	/*--------------------------------------------------------------------------------*/
	/* -------------------------- Parse OUTPUT  ------------------------------------- */
	/*--------------------------------------------------------------------------------*/
	/*--------------------------------------------------------------------------------*/
	
	
	
	numdimsx  = numdimsm;
	
	
	dimsx     = (int *)mxMalloc(numdimsx*sizeof(int));
	
	
	dimsx[0]  = N;
	
	
	for (i = 1 ; i < numdimsx ; i++)
		
	{
		
		dimsx[i] = dimsm[i];
		
	}
	
	/* ----- output 1 ----- */
	
	
	plhs[0]    = mxCreateNumericArray(numdimsx, dimsx, mxDOUBLE_CLASS, mxREAL);
	
	x          = mxGetPr(plhs[0]);
	
	
	/* ----- output 2 ----- */
	
	
	plhs[1]    = mxCreateNumericArray(numdimsx, dimsx, mxDOUBLE_CLASS, mxREAL);
	
	y          = mxGetPr(plhs[1]);
	
	
	/* ----- Usefull Matrices & vectors ----- */
	
	
	cos_angle        = (double *)mxMalloc(N*sizeof(double));
	
	sin_angle        = (double *)mxMalloc(N*sizeof(double));
	
	
	/*---------------------------------------------------------------------------------*/
	/*---------------------------------------------------------------------------------*/
	/* ----------------------- MAIN CALL  -------------------------------------------- */
	/*---------------------------------------------------------------------------------*/
	/*---------------------------------------------------------------------------------*/
	/*---------------------------------------------------------------------------------*/
	
	
	ellipse(m , e , N , 
		x , y,
		slice, 
		cos_angle , sin_angle);
	
	
	
	/*-----------------------------------------------*/
	/*-----------------------------------------------*/
	/* ------------ END of Mex File ---------------- */
	/*-----------------------------------------------*/
	/*-----------------------------------------------*/
	
	
	mxFree(dimsx);
	
	mxFree(cos_angle);
	
	mxFree(sin_angle);
	
	
}


/*----------------------------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------------------------*/


void ellipse(double *m , double *e , int N , 
			 double *x , double *y,
			 int slice, 
			 double *cos_angle , double *sin_angle)
			 
			 
{
	
	double twoPI = 2.0*PI , pas_PI;
	
	double mx , my , a , b ;
	
	
	double theta , costheta , sintheta ;
	
	double acostheta , bcostheta , asintheta , bsintheta;
	
	int i , v;
	
	int  v2 , v3 , vN , ivN;
	
	
	pas_PI = twoPI/(N - 1);
	
	for(i = 0 ; i < N ; i++)
		
	{
		
		cos_angle[i] = cos(i*pas_PI);
		
		sin_angle[i] = sin(i*pas_PI);
		
	}
	
	
	for (v = 0 ; v < slice ; v++)
		
	{
		
		
		v2          = v*2;
		
		v3          = v*3;
		
		vN          = v*N;
		
		
		mx          = m[0 + v2];
		
		my          = m[1 + v2];
		
		
		
        a           = e[0 + v3];
		
		b           = e[1 + v3];
		
		theta       = e[2 + v3];
		
		
		
		
		costheta    = cos(theta);
		
		sintheta    = sin(theta);
		
		acostheta   = a*costheta;
		
		asintheta   = a*sintheta;
		
		bcostheta   = b*costheta;
		
		bsintheta   = b*sintheta;
		
		
		for (i = 0 ; i < N ; i++)
			
		{
			
			ivN       = i + vN;
			
			
			x[ivN]    = mx + acostheta*cos_angle[i] - bsintheta*sin_angle[i]; 
			
			
			y[ivN]    = my + asintheta*cos_angle[i] + bcostheta*sin_angle[i]; 
			
		}
		
	}
	
	
}



