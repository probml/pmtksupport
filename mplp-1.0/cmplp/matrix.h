#include <iostream>
#include <fstream>

class Matrix;

class Vec
{
 public:
  int m_size;
  double *m_dat;
  double *m_ep;
  Vec(int size)
  {
    m_size = size;    
    if (size==0) return;
    m_dat= new double[size];
    m_ep = &m_dat[size];
  }
  ~Vec();
  inline int size() {return m_size;}
  double & operator[](int i) {return m_dat[i];}
  Vec & operator=(double val);
  Vec & operator=(Vec & v);
  Vec & operator*=(double val);
  Vec & operator/=(double val);	
  Vec & operator+=(double val);
  Vec & operator-=(double val);
  Vec & operator+=(Vec & v);
  Vec & operator*=(Vec & v);
  Vec & operator/=(Vec & v);
  void MatMult(Matrix & m, Vec & v);
  void OuterProd(Vec &v, Matrix &m);   // m should be size this->Size,v.Size
  Vec & operator-=(Vec & v);
  Vec & BackDiv(Vec &v);  
  double Max();
  double Min();
  void MaxMin(double &mx, double &mn);
  double Sum();
/*  Vec & Log();
  Vec & Exp();*/
/*  double KL_Dist(Vec & v);*/
  Vec & Inv();
  void SetSize(int n);
  void Rand();
  void Read(char *fname);
  void Write(char *fname);
  void Print();
};


class Matrix
{
 public:
  int m_nrows;
  int m_ncols;
  bool m_bDelete;
  double *m_raw_dat;
  double *m_ep;
  double **m_rows;
  Matrix(int rows =0,  int cols = 0);
  Matrix(int rows, int cols,double *p_dat);
  ~Matrix();
  void SetSize(int rows, int cols);
  inline int nrows() {return m_nrows;}
  inline int ncols() {return m_ncols;}
  inline double & operator()(int r,int c) {return m_rows[r][c];}
  Matrix & BackDiv(Matrix &m);  
  Matrix & operator+=(Matrix &m);
  Matrix & operator-=(Matrix &m);	
  Matrix & operator/=(Matrix &m);
  Matrix & operator*=(Matrix &m);
  void VecMult(Vec &m, Vec & v);
  Matrix & operator*=(double val);
  Matrix & operator/=(double val);
  Matrix & operator+=(double val);
  Matrix & operator=(double val);
  Matrix & operator=(Matrix &m);
  void Rand();
  void Read(char *fname);
  void Write(char *fname);
  void Exp();
/*  double KL_Dist(Matrix & p);
  void Log();*/
  double Sum();
  double *GetRow(int i) {return m_rows[i];}
  void GetCol(int i, Vec &v);
  void GetRow(int i, Vec &v);
  void SetRow(int i, Vec & v);
  double SumRow(int i);	
  double SumCol(int i);
  void SumRows(Vec &v);
  void SumCols(Vec & v);
  void ScaleRow(int i,double s);
  // Scale the i'th row by v[i]
  Matrix & ScaleRowsByVec(Vec &v);  
  Matrix & DivRowsByVec(Vec &v);  
  void Transpose( Matrix & dst);
  void Print();
};

class SparseMatrix
{
public:
	double *m_dat;
	int *m_rind;
	int *m_cind;
	int m_nrows;
	int m_ncols;
	int m_nelem;
	double m_thr;

	double GetElem(int i, int & row, int & col) {row=m_rind[i]; col=m_cind[i]; return m_dat[i];}
//	double KL_Dist(Matrix &p);
	void Read(char *fname, double thr);
	double SumRow(int i);
	double SumCol(int i);
};

inline double dot(double *a, double *b, int n)
{
  double s=0;
  double *ap,*bp;
  double *aend=&a[n];
  for (ap=&a[0],bp=&b[0]; ap<aend; ap++, bp++)
      s+=*ap*(*bp);
  return s;
}

void mult_trans1(Matrix & A, Matrix & B, Matrix & C);
void mult_trans2(Matrix & A, Matrix & B, Matrix & C);
void mult_trans3(Matrix & A, Matrix & B, Matrix & C);
void CalcModel(Matrix & phis, Matrix & psis, Vec & a, Vec &b, Matrix & Model);
void CalcModelTrans(Matrix & phis, Matrix & psis, Vec & a, Vec &b, Matrix & Model);
void MatMaxEnt( Matrix & f_val, Matrix & f_exp, Matrix & p_proj, int ep, Matrix & lambdas, Vec & zs);
void ConjGradMaxEnt( Matrix & f_val, Matrix & f_exp, Matrix & p_proj, int ep, Matrix & lambdas, Vec & zs);
void MakeCondDistX(Matrix & pygx, Vec & px);
void VecMaxEnt( Matrix & f_val, Matrix & f_exp, Vec & p0, int ep, Matrix & lambdas, Vec & zs, bool brows, double maxentthr);
double dot(double *a, double *b, int n);
void ConjGradMaxEntVec( Matrix & f_val, Matrix & ref_p, Matrix & p_proj, int ep, Matrix & lambdas, Vec & zs, bool brow);
