#include <stdlib.h>     
#include <math.h>      
#include "various.h"
#define PI 3.141592653589793
#define EE 2.718281828459046


double rexp(double lambda)
{return -log(drand48())/lambda;}

double rgamma(double alpha) /* rgamma(alpha)/beta=rbeta(alpha,beta) */
{
   double r1,r2,aa,x,w,c1,c2,c3,c4,c5;
   if(alpha<=0.)return 0.;
   if(alpha == 1.)return rexp(1.);
   if(alpha<1){
      aa=(alpha+EE)/EE;
      do{
         r1=drand48();
         r2=drand48();
         if(r1>1./aa){
            x = -log(aa*(1.-r1)/alpha);
            if(r2<pow(x,(alpha-1.)))return x;
         }
         else{
            x = pow((aa*r1),(1./alpha));
            if(r2<exp(-x))return x;
         }
      }while(r2<2);
   }
   else{
      c1=alpha-1;
      c2=(alpha-1./(6.*alpha))/c1;
      c3=2./c1;
      c4=c3+2.;
      c5=1./sqrt(alpha);
      do{
         do{
            r1=drand48();
            r2=drand48();
            if(alpha>2.5)r1=r2+c5*(1.-1.86*r1);
         }while(r1<=0 || r1 >= 1);
         w=c2*r2/r1;
         if(c3*r1+w+1/w <= c4)return c1*w;
         if(c3*log(r1)-log(w)+w<1)return c1*w;
      }while(r2<2);
   }
}
/*************************************/
double rchisq(double t)
{return rgamma(t/2)*2.;}
/*************************************/

/*************************************/
double rnor(double mu,double sd)
{
int ixxx1=0;
   double e,v1,v2,w, rxxx1;
   //if(ixxx1==0){
      do{
         v1=2*drand48()-1.;
         v2=2*drand48()-1.;
         w=v1*v1+v2*v2;      
      }while(w>1.);
      e=sqrt((-2.*log(w))/w);
      rxxx1=v1*e;
      ixxx1=1;
      return v2*e*sd+mu;
   //}
   //else{
   //   ixxx1=0;
   //   return rxxx1*sd+mu;
   //}
}
/*************************************/
