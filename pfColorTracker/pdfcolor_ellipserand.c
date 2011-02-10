
#include <math.h>
#include <time.h>
#include "mex.h"


/*


  Usage 
  -----


  [py , zi , yi]   = pdfcolor_ellipserand(Z , y , e , Npdf , vect_edge1 , vect_edge2 , vect_edge3 );%


  Inputs
  ------

   Z               Image (m x n x 3)
   y               State position (2 x N)
   e               Ellipsoid (3 x N)
   Npdf            Number of samples to copute the color histogram
   vect_edge1      Color support (R or H) (1 x Nx)
   vect_edge2      Color support (G or S) (1 x Ny)
   vect_edge3      Color support (B or V) (1 x Nz)


  Ouputs
  -------

  py               pdf color (NxNyNz x N)
  zi               Interpolated color values (3 x Npf x N)
  yi               Position of Interpolated values (2 x Npf x N)



  To compile
  -----------

  mex -output pdfcolor_ellipserand.dll pdfcolor_ellipserand.c

  mex -f mexopts_intel10.bat -output pdfcolor_ellipserand.dll pdfcolor_ellipserand.c

  
  Example 1
  ---------

	  Z                = rand(200 , 200 , 3);
	  y                = [111 , 52 , 15 , 56 ; 22 , 100 , 34 , 43];
	  e                = [10 , 10 , 10 , 30 ; 15, 3 , 22 , 10 ; -pi/1.1 , pi/2 , 0 ,0.75*pi];
	  Npdf             = 4000;
	  Nx               = 4;
	  Ny               = 3;
	  Nz               = 5;
	  M                = Nx*Ny*Nz;
      vect_edge1       = (0 : 1/Nx : 1);
      vect_edge2       = (0 : 1/Ny : 1);
      vect_edge3       = (0 : 1/Nz : 1);
      [py , zi , yi]   = pdfcolor_ellipserand(Z , y , e , Npdf , vect_edge1 , vect_edge2 , vect_edge3 );%
	  figure(1)
	  imagesc(Z);
	  hold on
	  plot(squeeze(yi(1 , : , :)) , squeeze(yi(2 , : , :)) , '+')
	  axis ij
	  axis equal
	  hold off
	  figure(2)
	  plot((1:M) , py);
	  axis([1 , M , 0 , 0.12])


   Author : Sébastien PARIS  © (sebastien.paris@lsis.org) (2004)
   -------

	  
*/


#define NUMERICS_FLOAT_MIN 1.0E-37
#define PI 3.14159265358979323846
#define NOBIN -1

#ifndef max
    #define max(a,b) (a >= b ? a : b)
    #define min(a,b) (a <= b ? a : b)
#endif


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
#define randint SHR3
#define rand() (0.5 + (signed)randint*2.328306e-10)



#ifdef __x86_64__
    typedef int UL;
#else
    typedef unsigned long UL;
#endif


/*--------------------------------------------------------------------------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------------------------------*/


static UL jsrseed = 31340134 , jsr;



/*-------------------------------------------------------------------------*/
/*-------------------------------------------------------------------------*/
/*------------------------ Headers   --------------------------------------*/
/*-------------------------------------------------------------------------*/
/*-------------------------------------------------------------------------*/

void  randini(void);  
int bin(double , double * , int);

/*-------------------------------------------------------------------------*/
/*-------------------------------------------------------------------------*/



void mexFunction( int nlhs, mxArray *plhs[] , int nrhs, const mxArray *prhs[] )

{
	
	double *Z , *y , *e  ;
	
	double *py , *zi , *yi;
	
	double *vect_edge1 , *vect_edge2 , *vect_edge3;
	
	double cte1 , cte2 , cte3 , twoPI = 2*PI , u , angle , val ;
	
	double theta , cos_theta , sin_theta , tmp1 , tmp2 , s , t , invM;

	//double NAN; NAN = mxGetNaN(); 
    //double NAN = mxGetNaN();
	
	
	const int *dimsZ=NULL , *dimsy=NULL , *dimse=NULL , *dimsvect_edge1=NULL , *dimsvect_edge2=NULL , *dimsvect_edge3=NULL;
	
	int *dimspy=NULL , *dimszi=NULL , *dimsyi=NULL;
	
	int i , j , n , m , nm , r , N , Npdf , NNpdf , Nx , Ny , Nz , Nxy;
	
	int numdimsZ , numdimsy , numdimse , numdimsvect_edge1 , numdimsvect_edge2 , numdimsvect_edge3;
	
	int numdimspy , numdimszi , numdimsyi; 
	
	int iNpdf2, iNpdf3 , iM , i2 , i3 , j2 , j3 , ndx , kx , ky , kz , M , j2iNpdf2 , j3iNpdf3;
	
	int floors , floort , compteur;


	
	
	/*---------------------------------------------------------------*/
	/*---------------------- PARSE INPUT ----------------------------*/	
	/*---------------------------------------------------------------*/
	
	if(nrhs < 3)
		
	{
		mexErrMsgTxt("At least 3 inputs are required");		
	}
	
	/*----------- Input 1 ---------------*/
	
    
    Z        = mxGetPr(prhs[0]);
    
    numdimsZ = mxGetNumberOfDimensions(prhs[0]);
    
	dimsZ    = mxGetDimensions(prhs[0]);
	
    if (numdimsZ != 3)
		
	{
		mexErrMsgTxt("First input must at least a Matrix (n x m x 3)");
	}
	
	
	n  = dimsZ[0];
	
	m  = dimsZ[1];
	
	nm = n*m;
	
	
	/*----------- Input 2 ---------------*/
	
	y         = mxGetPr(prhs[1]);
    
    numdimsy  = mxGetNumberOfDimensions(prhs[1]);
    
	dimsy     = mxGetDimensions(prhs[1]);
	
	if (numdimsy != 2)
		
	{	
		mexErrMsgTxt("Second input be (2 x N)");
	}
	
	if (dimsy[0] != 2)
		
	{
		
		mexErrMsgTxt("Second input be (2 x N)");
		
	}
	
	N  = dimsy[1];
	
	/*----------- Input 3 ---------------*/
	
	
	e        = mxGetPr(prhs[2]);
    
    numdimse = mxGetNumberOfDimensions(prhs[2]);
    
	dimse    = mxGetDimensions(prhs[2]);
	
	if (numdimse != 2)
		
	{	
		mexErrMsgTxt("Third input be (3 x N)");
	}
	
	if (dimse[0] != 3 || dimse[1] != N)
		
	{
		
		mexErrMsgTxt("Second input be (3 x N)");
		
	}

	/*----------- Input 4 ---------------*/

		if(nrhs >= 4)
		
	{
		
		Npdf = (int)mxGetScalar(prhs[3]);
		
	}
	
	else
		
	{
		
		Npdf = 100;
	}
	
	
	NNpdf = Npdf*N;

	/*----------- Input 5 ---------------*/

	
	if(nrhs >= 5)
		
	{
		
		vect_edge1        = mxGetPr(prhs[4]);
		
		numdimsvect_edge1 = mxGetNumberOfDimensions(prhs[4]);
		
		dimsvect_edge1    = mxGetDimensions(prhs[4]);
		
		if (numdimsvect_edge1 != 2)
			
		{
			
			mexErrMsgTxt("Fifth input be vector");	
			
		}
		
		if (((dimsvect_edge1[0] != 1) && (dimsvect_edge1[1] > 1)) || ((dimsvect_edge1[1] != 1) && (dimsvect_edge1[0] > 1)))
			
		{	
			mexErrMsgTxt("Fourth input be (1 x (Nx + 1)) or ((Nx + 1)x1)");	
		}
		
        Nx    = max(dimsvect_edge1[0] , dimsvect_edge1[1]) - 1;
		
	}
	
	else
		
	{
		
		vect_edge1 = (double *)mxMalloc(1*9*sizeof(double));
		
        Nx         = 8;
		
		cte1       = 1/Nx;
		
		for (i = 0; i < Nx + 1 ; i++)
			
		{		
			vect_edge1[i] = i*cte1;
		}
				
	}
	
	/*----------- Input 6 ---------------*/

	
	if(nrhs >= 6)
		
	{
		
		vect_edge2        = mxGetPr(prhs[5]);
		
		numdimsvect_edge2 = mxGetNumberOfDimensions(prhs[5]);
		
		dimsvect_edge2    = mxGetDimensions(prhs[5]);
		
		if (numdimsvect_edge2 != 2)
			
		{
			
			mexErrMsgTxt("Sixth input be vector");	
			
		}
		
		if (((dimsvect_edge2[0] != 1) && (dimsvect_edge2[1] > 1)) || ((dimsvect_edge2[1] != 1) && (dimsvect_edge2[0] > 1)))
			
		{	
			mexErrMsgTxt("Sixth input be (1 x (Nx + 1)) or ((Nx + 1)x1)");	
		}
		
        Ny    = max(dimsvect_edge2[0] , dimsvect_edge2[1]) - 1;
		
	}
	
	else
		
	{
		
		vect_edge2 = (double *)mxMalloc(1*9*sizeof(double));
		
        Ny         = 8;
		
		cte2       = 1/Ny;
		
		for (i = 0; i < Ny + 1 ; i++)
			
		{		
			vect_edge2[i] = i*cte2;
		}
				
	}

	/*----------- Input 7 ---------------*/

	
	if(nrhs >= 7)
		
	{
		
		vect_edge3        = mxGetPr(prhs[6]);
		
		numdimsvect_edge3 = mxGetNumberOfDimensions(prhs[6]);
		
		dimsvect_edge3    = mxGetDimensions(prhs[6]);
		
		if (numdimsvect_edge3 != 2)
			
		{
			
			mexErrMsgTxt("Seven input be vector");	
			
		}
		
		if (((dimsvect_edge3[0] != 1) && (dimsvect_edge3[1] > 1)) || ((dimsvect_edge3[1] != 1) && (dimsvect_edge3[0] > 1)))
			
		{	
			mexErrMsgTxt("Seven input be (1 x (Nx + 1)) or ((Nx + 1)x1)");	
		}
		
        Nz    = max(dimsvect_edge3[0] , dimsvect_edge3[1]) - 1;
		
	}
	
	else
		
	{
		
		vect_edge3 = (double *)mxMalloc(1*9*sizeof(double));
		
        Nz         = 8;
		
		cte3       = 1/Nz;
		
		for (i = 0; i < Nz + 1 ; i++)
			
		{		
			vect_edge3[i] = i*cte3;
		}
				
	}


	M        = Nx*Ny*Nz;
	
	Nxy      = Nx*Ny;
	
	invM     = 1.0/(double)(M);
	
	
	
	/*----------- Output 1 ---------------*/
	
	numdimspy      = 2;
	
	dimspy         = (int *)mxMalloc(numdimspy*sizeof(int));
	
	dimspy[0]      = M;
	
    dimspy[1]      = N;
	
	
	
	plhs[0]        = mxCreateNumericArray(numdimspy, dimspy, mxDOUBLE_CLASS, mxREAL);
	
	py             = mxGetPr(plhs[0]);
	
	
	/*----------- Output 2 ---------------*/
	
	
	numdimszi      = 3;
	
	dimszi         = (int *)mxMalloc(numdimszi*sizeof(int));
	
	dimszi[0]      = 3;
	
    dimszi[1]      = Npdf;
	
    dimszi[2]      = N;
	
	
	
	plhs[1]        = mxCreateNumericArray(numdimszi, dimszi, mxDOUBLE_CLASS, mxREAL);
	
	zi             = mxGetPr(plhs[1]);
	
	
	/*----------- Output 3 ---------------*/
	
	
	numdimsyi      = 3;
	
	dimsyi         = (int *)mxMalloc(numdimsyi*sizeof(int));
	
	dimsyi[0]      = 2;
	
    dimsyi[1]      = Npdf;
	
    dimsyi[2]      = N;
	
	
	
	plhs[2]        = mxCreateNumericArray(numdimsyi, dimsyi, mxDOUBLE_CLASS, mxREAL);
	
	yi             = mxGetPr(plhs[2]);

	
	
	
	/* 1) --------- Compute ellipse Coordinates & Interpolate values & Histograms---------------*/
	
	
	randini();


	for (i = 0 ; i < N ; i++)
		
	{
		
		iNpdf2     = i*Npdf*2;
		
		iNpdf3     = i*Npdf*3;
		
		iM         = i*M;
		
		i2         = i*2;
		
		i3         = i*3;
		
		theta      = e[2 + i3];
		
		cos_theta  = cos(theta);
		
		sin_theta  = sin(theta);
		
		compteur = 0;
		
		
		for (j = 0 ; j < Npdf ; j++)
			
		{
			
			j2                  = 2*j;

			j2iNpdf2            = j2 + iNpdf2;
			
    	    j3                  = 3*j;

			j3iNpdf3            = j3 + iNpdf3;
			
			u                   = rand();

			angle               = twoPI*rand();

			
			tmp1                = (e[0 + i3]*u)*cos(angle);
			
			tmp2                = (e[1 + i3]*u)*sin(angle);
			
			s                   = y[0 + i2] + (cos_theta*tmp1 - sin_theta*tmp2);
			
			t                   = y[1 + i2] + (sin_theta*tmp1 + cos_theta*tmp2);
			
			yi[0 + j2iNpdf2]    =  s;
			
			yi[1 + j2iNpdf2]    =  t;
			
			
			if ((s < 1) || (s > m) || (t < 1) || (t > n))
				
			{
				

				zi[0 + j3iNpdf3]    = NAN;
				
				zi[1 + j3iNpdf3]    = NAN;
				
				zi[2 + j3iNpdf3]    = NAN;
				
				kx                  = NOBIN;
				
				ky                  = NOBIN;
				
				kz                  = NOBIN;
				
			}
			
			else
				
			{
				
                floors              = floor(s);
				
				floort              = floor(t);
				
				ndx                 = floort + (floors - 1)*n - 1;
				
				
				s                   = (s - floors);
				
				t                   = (t - floort);
				
				
				r                   = 0*nm;
				
				val                 = (Z[ndx + r]*(1 - t) + Z[ndx + r + 1]*t)*(1 - s) + ( Z[ndx + r + n]*(1 - t) + Z[ndx + n + 1 + r]*t )*s;
				
				zi[0 + j3iNpdf3]    = val;
				
                kx                  =  bin(val , vect_edge1 , Nx);
				
				
				
				r                   = 1*nm;
				
				val                 = (Z[ndx + r]*(1 - t) + Z[ndx + r + 1]*t)*(1 - s) + ( Z[ndx + r + n]*(1 - t) + Z[ndx + n + 1 + r]*t )*s;
				
				zi[1 + j3iNpdf3]    = val;
				
				ky                  =  bin(val , vect_edge2 , Ny);
				
				
				
				r                   = 2*nm;
				
				val                 = (Z[ndx + r]*(1 - t) + Z[ndx + r + 1]*t)*(1 - s) + ( Z[ndx + r + n]*(1 - t) + Z[ndx + n + 1 + r]*t )*s;
				
				zi[2 + j3iNpdf3]    = val;
				
				kz                  =  bin(val , vect_edge3 , Nz);


				
			}
			
			if ((kx != NOBIN) && (ky != NOBIN) && (kz != NOBIN))
				
			{
				
				py[kx + ky*Nx + kz*Nxy + iM]++;
				
				compteur++;
				
			}
			
			
		}
		
		if (compteur != 0)
		{
			
			invM        = 1.0/((double)compteur);
			
		}
		
		else
			
		{
			
			invM = 1.0;
			
		}
		for (j = 0 ; j < M ; j++)
			
		{
						
			py[j + iM] *= invM;
			
		}
		
		
	}
	
	
	
	mxFree(dimspy);
	
	mxFree(dimszi);
	
	mxFree(dimsyi);
	
		
	
}


/*------------------------------------------------------------*/

int bin(double zi , double *vect_edge , int Nbin)

{
	
	int k  = NOBIN;
	
	int k0 = 0 , k1 , N1 = Nbin;
	
	k1     = N1;
	
	
	if ((zi >= vect_edge[0]) && (zi <= vect_edge[N1]))
	{
		
		k = (k0 + k1)/2;
		
		while (k0 < k1 - 1)
			
		{
			if (zi >= vect_edge[k]) 
				
			{
				k0 = k;
			}
			
			else 
			{
				k1 = k;
			}
			
			k = (k0 + k1)/2;
		}
		
		k = k0;
	}
	
	return k;
	   
}

 /* ----------------------------------------------------------------------- */

void randini(void) 

{
	
	/* SHR3 Seed initialization */
	
	jsrseed  = (UL) time( NULL );
	
	jsr     ^= jsrseed;
	
		
	
}


