#include "matrix.h"

void CheckPomwSizeer(void* pomwSizeer);
mwSize choose(mwSize m, mwSize n);
double logdexp(double x,double beta);
double logdpois(double x,double lambda);
double getLogBetaPrior(double m,double x);
double** ReadData();
void FreeData(double** Data);
void NormalizeData(double** data);
double getSSD(double** data,mwSize vi,mwSize vj);
double mylgamma(double lambda,mwSize p);
bool FileExists(char* strFilename);
#if defined(_WIN32)
double lgamma(double x);
double drand48();
#endif
#if defined(_WIN64)
double lgamma(double x);
double drand48();
#endif
