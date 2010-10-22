#include "matrix.h"

using namespace std;

Vec::~Vec()
{
  delete [] m_dat;
}

Vec & Vec::operator=(double val)
{
  for (double *p=m_dat ; p<m_ep; p++)
	  (*p)=val;
  return (*this);
}

Vec & Vec::operator=(Vec & v)
{
  double *pv = v.m_dat;
  for (double *p=m_dat ; p<m_ep; p++,pv++)
	  (*p)=(*pv);
  return (*this);
}

void Vec::Read(char *fname)
{
	ifstream ifs(fname);
	int s;

	ifs >> s;
	if (m_dat)
	   delete [] m_dat;
	m_size = s;
	m_dat= new double[m_size];
	m_ep = &m_dat[m_size];
    for (int i=0; i< s; i++)
		ifs >> m_dat[i];
}

void Vec::Write(char *fname)
{
	ofstream ofs(fname);

	ofs << size() << endl;
    for (int i=0; i< size(); i++)
		ofs << m_dat[i] << " ";
}

void Vec::SetSize(int n)
{
	if (m_dat)
	   delete [] m_dat;
	m_size = n;
	m_dat= new double[m_size];
	m_ep = &m_dat[m_size];
}

void Vec::Rand()
{
  for (int i=0; i< size(); i++)
      m_dat[i]=(double) 2.0*rand()/RAND_MAX -1;
}

double Vec::Max()
{
  double m = m_dat[0];
  for (int i=1;i<size(); i++) {
    if (m_dat[i]>m)
      m=m_dat[i];
  }
  return m;
}

double Vec::Min()
{
  double m = m_dat[0];
  for (int i=1;i<size(); i++) {
    if (m_dat[i]<m)
      m=m_dat[i];
  }
  return m;
}

void Vec::MaxMin(double &mx, double &mn)
{
	mx = m_dat[0];
	mn = mx;
	for (int i=1;i<size(); i++) 
	{
		if (m_dat[i]>mx)
			mx=m_dat[i];
		if (m_dat[i]<mn)
			mn=m_dat[i];
	}
}


double Vec::Sum()
{
  double s=0;
   for (double *p=m_dat; p<m_ep; p++)
	   s+=(*p);
  return s;
}

void Vec::Print()
{
  for (int i=0;i<size(); i++) 
	  cout << m_dat[i] << " ";
  cout << endl;
}

/*
Vec & Vec::Log()
{
   for (double *p=m_dat; p<m_ep; p++)
	   (*p)=(double) log((*p));
  return (*this);
}

Vec & Vec::Exp()
{
   for (double *p=m_dat; p<m_ep; p++)
	   (*p)=(double) exp((*p));
  return (*this);
}
*/

Vec & Vec::operator*=(double val)
{
   for (double *p=m_dat; p<m_ep; p++)
	   (*p)*=val;
  return (*this);
}

Vec & Vec::operator/=(double val)
{
   for (double *p=m_dat; p<m_ep; p++)
	   (*p)/=val;
  return (*this);
}

Vec & Vec::Inv()
{
   for (double *p=m_dat; p<m_ep; p++)
	   (*p)=1/(*p);
  return (*this);
}

Vec & Vec::operator+=(double val)
{
   for (double *p=m_dat; p<m_ep; p++)
	   (*p)+=val;
  return (*this);
}

Vec & Vec::operator-=(double val)
{
   for (double *p=m_dat; p<m_ep; p++)
	   (*p)-=val;

  return (*this);
}

Vec & Vec::operator+=(Vec & v)
{
	double *p,*pv;
   for (p=m_dat, pv=v.m_dat; p<m_ep; p++, pv++)
	   (*p)+=(*pv);
  return (*this);
}

Vec & Vec::operator-=(Vec & v)
{
   for (double *p=m_dat, *pv=v.m_dat; p<m_ep; p++, pv++)
	   (*p)-=(*pv);
  return (*this);
}

Vec & Vec::operator*=(Vec & v)
{
   for (double *p=m_dat, *pv=v.m_dat; p<m_ep; p++, pv++)
	   (*p)*=(*pv);
  return (*this);
}

Vec & Vec::operator/=(Vec & v)
{
   for (double *p=m_dat, *pv=v.m_dat; p<m_ep; p++, pv++)
	   (*p)/=(*pv);
  return (*this);
}

Vec & Vec::BackDiv(Vec &v)
{
   for (double *p=m_dat, *pv=v.m_dat; p<m_ep; p++, pv++)
	   (*p)=(*pv)/(*p);

  return (*this);

}

void Vec::MatMult(Matrix & m, Vec & v)
{
  double *pv;
  double **mrows = m.m_rows;
  for (int i=0; i<m.m_ncols; i++)
    {
      double sum=0;
      pv = m_dat;
      for (int j=0; j<m.m_nrows; j++,pv++)
	sum+=(*pv)*mrows[j][i];	
	//	sum+=(*pv)*m(j,i);
      v[i]=sum;
    }
}

/*
double Vec::KL_Dist(Vec & v)
{
  double d=0;
  for (double *p=m_dat, *pv=v.m_dat; p<m_ep; p++, pv++)
    {
      if (!(*pv))
	*pv=(double) 1e-70;
      d+=(*p)* (double) log((*p)/(*pv));
    }
  return d;
}
*/

Matrix::Matrix(int rows, int cols) : m_nrows(rows), m_ncols(cols),m_bDelete(true)
{
  m_raw_dat = new double[rows*cols];
  m_rows = new double*[rows];
  m_ep = &m_raw_dat[rows*cols];
  for (int i=0; i<rows; i++)
    m_rows[i]=&m_raw_dat[i*cols];
}

Matrix::~Matrix()
{
  if (m_bDelete)
	  delete [] m_raw_dat;
  delete [] m_rows;

}

Matrix::Matrix(int rows, int cols,double *p_dat)  :  m_nrows(rows), m_ncols(cols),m_bDelete(false)
{
    m_raw_dat=p_dat;
	m_rows = new double*[rows];
    m_ep = &m_raw_dat[rows*cols];
    for (int i=0; i<rows; i++)
         m_rows[i]=&m_raw_dat[i*cols];

}

void Matrix::SetSize(int rows, int cols)
{
  this->~Matrix();
  
  m_nrows=rows;
  m_ncols=cols;


  m_raw_dat = new double[rows*cols];
  m_ep = &m_raw_dat[rows*cols];
  m_rows = new double*[rows];
  
  for (int i=0; i<rows; i++)
    m_rows[i]=&m_raw_dat[i*cols];
}

void Matrix::SetRow(int i, Vec & v)
{
	memcpy(m_rows[i],v.m_dat,m_ncols*sizeof(double));
}

void Matrix::GetRow(int i, Vec &v)
{
	memcpy(v.m_dat,m_rows[i],m_ncols*sizeof(double));
}

/*
double Matrix::KL_Dist(Matrix & m)
{
  double d=0;
   for (double *p=m_raw_dat, *pm=m.m_raw_dat; p<m_ep; p++, pm++)
	   d+=(*p)* (double) log((*p)/(*pm));

  return d;
}
*/

void Matrix::Rand()
{
  for (int i=0; i< m_nrows; i++)
    for (int j=0; j< m_ncols; j++)
      (*this)(i,j)=(double) 2.0* (double) rand()/RAND_MAX -1;
}

void Matrix::Read(char *fname)
{
	ifstream ifs(fname);
	int rows,cols;

	ifs >> rows;
	ifs >> cols;
	SetSize(rows,cols);
	for (int i=0; i< m_nrows; i++)
	  for (int j=0; j< m_ncols; j++)
	    ifs >> (*this)(i,j);
}


void Matrix::Write(char *fname)
{
	ofstream ofs(fname);

	ofs << m_nrows << " ";
	ofs << m_ncols << " " << endl;
    for (int i=0; i< m_nrows; i++)
	{
       for (int j=0; j< m_ncols; j++)
		   ofs << (*this)(i,j) << " ";
	   ofs << endl;
	}
}

Matrix & Matrix::operator+=(Matrix &m)
{
  for (double *p=m_raw_dat,*pm=m.m_raw_dat; p<m_ep; p++,pm++)  
	  (*p)+=(*pm);
  return (*this);
}

Matrix & Matrix::operator-=(Matrix &m)
{
  for (double *p=m_raw_dat,*pm=m.m_raw_dat; p<m_ep; p++,pm++)  
	  (*p)-=(*pm);
  return (*this);
}

Matrix & Matrix::operator/=(Matrix &m)
{
  for (double *p=m_raw_dat,*pm=m.m_raw_dat; p<m_ep; p++,pm++)  
	  (*p)/=(*pm);
  return (*this);
}

Matrix & Matrix::operator*=(Matrix &m)
{
  for (double *p=m_raw_dat,*pm=m.m_raw_dat; p<m_ep; p++,pm++)  
	  (*p)*=(*pm);
  return (*this);
}

Matrix & Matrix::operator+=(double val)
{
  for (double *p=m_raw_dat; p<m_ep; p++)
	  (*p)+=val;
  return (*this);
}

Matrix & Matrix::BackDiv(Matrix &m)
{
  for (int i=0; i< m_nrows; i++)
    for (int j=0; j< m_ncols; j++)
      (*this)(i,j)=m(i,j)/(*this)(i,j);
  return (*this);
}

void Matrix::VecMult(Vec &m, Vec & v)
{
  double *pv=v.m_dat;
  double **pm=m_rows;
  
  for (int i=0; i<m_nrows; i++,pv++,pm++)
    (*pv)=dot(m.m_dat,(*pm),m_ncols);
}

Matrix & Matrix::operator*=(double val)
{
  for (double *p=m_raw_dat; p<m_ep; p++)
	  (*p)*=val;
  return (*this);
}

Matrix & Matrix::operator/=(double val)
{
  for (double *p=m_raw_dat; p<m_ep; p++)
	  (*p)/=val;
  return (*this);
}

Matrix & Matrix::operator=(double val)
{
  for (double *p=m_raw_dat; p<m_ep; p++)
	  (*p)=val;
  return (*this);
}

Matrix & Matrix::operator=(Matrix &m)
{
  for (double *p=m_raw_dat,*pm=m.m_raw_dat; p<m_ep; p++,pm++)  
	  (*p)=(*pm);
  return (*this);
}


/*
void Matrix::Exp()
{
  for (double *p=m_raw_dat; p<m_ep; p++)
	  (*p)=(double) exp(*p);
}

void Matrix::Log()
{
  for (double *p=m_raw_dat; p<m_ep; p++)
	  (*p)=(double) log(*p);
}
*/

double Matrix::Sum()
{
  double s=0;
  for (double *p=m_raw_dat; p<m_ep; p++)
	  s+=(*p);
  return s;
}

double Matrix::SumRow(int i)
{
  double s=0;
  double *ep = &(m_rows[i][m_ncols]);

  ofstream ofs("c:\\dsdr\\ArData\\dbgme.txt");
  for (double *p=m_rows[i]; p<ep; p++)
	  s+=(*p);
  return s;
}

double Matrix::SumCol(int i)
{
	double sum=0;
	for (int j=0; j<m_nrows; j++)
		sum+=(*this)(j,i);
	return sum;
}

void Matrix::GetCol(int i, Vec &v)
{
	for (int j=0; j<m_nrows; j++)
		v[j]=(*this)(j,i);
}

void Matrix::ScaleRow(int i,double s)
{
  double *ep = &(m_rows[i][m_ncols]);

  for (double *p=m_rows[i]; p<ep; p++)
	  (*p)*=s;
}

Matrix & Matrix::ScaleRowsByVec(Vec &v)
{
  for (int i=0; i< m_nrows; i++)
    ScaleRow(i,v[i]);
  return (*this);
}

Matrix & Matrix::DivRowsByVec(Vec &v)
{
  for (int i=0; i< m_nrows; i++)
    ScaleRow(i,1/v[i]);
  return (*this);
}

void Matrix::SumRows(Vec &v)
{
	for (int i=0; i<m_nrows; i++)
		v[i]=SumRow(i);
}

void Matrix::SumCols(Vec & v)
{
	for (int i=0; i<m_ncols; i++)
		v[i]=SumCol(i);
}

void Matrix::Print()
{
  for (int i=0; i< m_nrows; i++)
    {
      for (int j=0; j< m_ncols; j++)
	cout << (*this)(i,j) << " ";
      cout << endl;
    }
}

void Matrix::Transpose( Matrix & dst)
{
  for (int i=0; i<m_nrows; i++)
    for (int j=0; j<m_ncols; j++)
      dst(j,i)=(*this)(i,j);
}

/*
double SparseMatrix::KL_Dist(Matrix &p)
{
	double sum=0;

	for (int i=0; i<m_nelem; i++)
	{
		double v1=m_dat[i];
		double v2=p(m_rind[i],m_cind[i]);
		sum+=v1* (double) log(v1/v2);
	}
	return sum;
}
*/

void SparseMatrix::Read(char *fname, double thr)
{
	ifstream ifs(fname);
	int rows,cols,nelem;

	ifs >> rows;
	ifs >> cols;
	m_nrows=rows;
	m_ncols=cols;
	ifs >> nelem;
	m_nelem=nelem;
	m_dat = new double[nelem];
	m_rind = new int[nelem];
	m_cind = new int[nelem];
	m_thr=thr;

	cout << thr << endl;
	double Z=0;
	int i;

    for (i=0; i< nelem; i++)
	{
		int row,col;
		double val;
		ifs >> row;
		ifs >> col;
		ifs >> val;
		m_dat[i]=val;
		m_rind[i]=row-1;
		m_cind[i]=col-1;
		Z+=val;
	}
	Z+=m_thr*(rows*cols-nelem);
	for (i=0; i<m_nelem; i++)
		m_dat[i]/=Z;
	m_thr/=Z;
}

double SparseMatrix::SumCol(int i)
{
	double sum=0;
	int count=0;
	for (int el=0; el<m_nelem; el++)
		if (m_cind[el]==i)
		{
			sum+=m_dat[el];
			count++;
		}
	sum+=(m_nrows-count)*m_thr;
	return sum;
}

double SparseMatrix::SumRow(int i)
{
	double sum=0;
	int count=0;
	for (int el=0; el<m_nelem; el++)
		if (m_rind[el]==i)
		{
			sum+=m_dat[el];
			count++;
		}
	sum+=(m_ncols-count)*m_thr;
	return sum;
}
