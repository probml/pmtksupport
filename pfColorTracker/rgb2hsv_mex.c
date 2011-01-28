#include <math.h>
#include "mex.h"

/*


  Usage: A    = rgb2hsv_mex(im)
  -----

  Inputs
  ------
 
  im         Image (n x m x 3) double

	
  Ouputs
  -------
	  
  Z          HSV image (n x m x 3)



  Example 
  -------


   im   = double(ceil(255*rand(200 , 300 , 3)));

   Z    = rgb2hsv_mex(im);


  To compile
  -----------

  mex rgb2hsv_mex.c

  mex -f mexopts_intel10.bat -output rgb2hsv_mex.dll rgb2hsv_mex.c



  Author : Sébastien PARIS  © (sebastien.paris@lsis.org) (2004)
  -------

  
*/

#define MIN(a , b) ((a) < (b) ? (a) : (b))
#define MAX(a , b) ((a) >= (b) ? (a) : (b))
#define NO_HUE   0


/*-------------------------------------------------------------------------*/
/*-------------------------------------------------------------------------*/
/*------------------------ Headers   --------------------------------------*/
/*-------------------------------------------------------------------------*/
/*-------------------------------------------------------------------------*/


void rgb2hsv (double * , double * , int  , int);


/*-------------------------------------------------------------------------*/
/*-------------------------------------------------------------------------*/



void mexFunction( int nlhs, mxArray *plhs[] , int nrhs, const mxArray *prhs[] )

{
	

	
	double *RGB; 
		
	double *HSV;
	
	const int *dimsRGB;
	
	int  n , m , numdimsRGB;
	
	
	/*---------------------------------------------------------------*/
	/*---------------------- PARSE INPUT ----------------------------*/	
	/*---------------------------------------------------------------*/
	
	if(nrhs != 1)
		
	{
		
		mexErrMsgTxt("1 input is required");
		
	}
	
	/*---------------- Input 1 ----------------------*/
	
    
    RGB        = mxGetPr(prhs[0]);
    
    numdimsRGB = mxGetNumberOfDimensions(prhs[0]);
    
	dimsRGB    = mxGetDimensions(prhs[0]);
	
    if (numdimsRGB != 3)
		
	{
		mexErrMsgTxt("First input must at least a Matrix (n x m x 3)");
	}
	
	
	n       = dimsRGB[0];
	
	m       = dimsRGB[1];
	
	
	/*------------------- Output 1 ----------------------*/
	
	
	plhs[0] = mxCreateNumericArray(numdimsRGB, dimsRGB, mxDOUBLE_CLASS, mxREAL);
	
	HSV     = mxGetPr(plhs[0]);
	
	
	/*----------- Convert into HSV ---------------*/
	
	
	rgb2hsv(RGB , HSV , n , m);
	
	
	/*--------------------------------------------*/
	
}


/* ----------------------------------------------------------------------------------------------- */
/* ----------------------------------------------------------------------------------------------- */



void rgb2hsv (double *RGB , double *HSV , int n ,  int m)
{
	
	
    double  r , g , b , max , min , delta;
	
	double  h , s , v , cte = 1.0/255.0, cte1 = 1.0/360.0;
	
	int     nm = n*m , i;
	
	
	for (i = 0; i < nm ; i++)
		
	{
		
		r     = RGB[i]*cte;
		
		g     = RGB[i + nm]*cte;
		
		b     = RGB[i + 2*nm]*cte;
		
		
		max   = MAX(r , MAX(g , b));
		
		min   = MIN(r , MIN(g , b));
		
		delta = (max - min);
		
		v     = max;
		
		
		if (max != 0.0)
			
		{
			s = delta/max;
		}
		
		else
			
		{
			s = 0.0;
		}
		
		if (s != 0.0) 
			
		{
			if (r == max)
				
			{
				h = (g - b)/delta;
			}
			
			else if (g == max)
				
			{
				h = 2.0 + (b - r)/delta;
			}
			
			else if (b == max)
				
			{
				h = 4.0 + (r - g)/delta;
			}
			
			h *= 60.0;
			
			if (h < 0) 
				
			{
				h += 360.0;
			}
			
			h *= cte1;			
			
		}
		
		else 
			
		{
			
			h = NO_HUE;
			
		}
		
		HSV[i]        = h;
		
		HSV[i + nm]   = s;
		
		HSV[i + 2*nm] = v;
		
	}
}

/* ----------------------------------------------------------------------------------------------- */
/* ----------------------------------------------------------------------------------------------- */
