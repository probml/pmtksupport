/*
 *  muldim_arr.cpp
 *  mplp
 *
 *  Created by Amir Globerson on 6/24/08.
 *  Copyright 2008 MIT. All rights reserved.
 *
 */
using namespace std;
#include "muldim_arr.h"

extern void print_int_vec(vector<int> v);

MulDimArr::MulDimArr(vector<int> & base_sizes)
{
	m_base_sizes = base_sizes;
	m_n_prodsize = 1;
	for (int i=0; i<base_sizes.size(); i++)
		m_n_prodsize*= base_sizes[i];
    m_dat= new double[m_n_prodsize];

	// Initialize array to zero
	for(int i=0; i< m_n_prodsize; i++)
		m_dat[i] = 0.0;
	
    m_ep = &m_dat[m_n_prodsize];	
};

MulDimArr::MulDimArr(const MulDimArr & v)
{
	m_base_sizes = v.m_base_sizes;
	m_n_prodsize = v.m_n_prodsize;
    m_dat= new double[m_n_prodsize];
    m_ep = &m_dat[m_n_prodsize];	
	memcpy(m_dat,v.m_dat,m_n_prodsize*sizeof(double));
}

void MulDimArr::print()
{
	for (int i=0; i<m_n_prodsize; i++)
	{
		cout << m_dat[i] << " ";
	}
	cout << endl;
}

// Print elements along with their indices
void MulDimArr::print_with_inds()
{
	vector<int> inds_for_big;
	int i;
	
	for (i=0; i< m_base_sizes.size(); i++)
		inds_for_big.push_back(0);
			
	for (i=0; i<m_n_prodsize; i++)
	{
		cout << "[";
		print_int_vec(inds_for_big);
		cout << "]: ";
		cout << m_dat[i] << endl;
		BaseInc(inds_for_big);
	}
}

MulDimArr & MulDimArr::operator=(const MulDimArr & v)
{
	m_base_sizes = v.m_base_sizes;
	m_n_prodsize = v.m_n_prodsize;
	if (m_dat!=NULL)
		delete [] m_dat; 
    m_dat= new double[m_n_prodsize];
    m_ep = &m_dat[m_n_prodsize];	
	memcpy(m_dat,v.m_dat,m_n_prodsize*sizeof(double));
	return (*this);
}

// Return the index in the flat multi-dimensional array corresponding to
// the given multi-index
inline int MulDimArr::GetFlatInd(vector<int> base_inds)
{
	int y=0;
	int fact = 1;
	int i;
	int nx = m_base_sizes.size();
	
	for (i=nx-1;i>=0; i--) {
		y += base_inds[i]*fact;
		fact*= m_base_sizes[i];
	}
	return y;
}


inline int MulDimArr::GetFlatIndFromBig(vector<int> big_base_inds, vector<int> inds_in_big)
{
	int y=0;
	int fact = 1;
	int i;
	int nx = m_base_sizes.size();
	
	for (i=nx-1;i>=0; i--) {
		y += big_base_inds[inds_in_big[i]]*fact;
		fact*= m_base_sizes[i];
	}
	return y;
}



double MulDimArr::GetVal(vector<int> & indices)
{
	return m_dat[GetFlatInd(indices)];
}


inline void MulDimArr::BaseInc(vector<int> & inds)
{
	int nx = m_base_sizes.size();
	int pos = nx-1;
	inds[pos]++;
	while (inds[pos]==m_base_sizes[pos]) {
		inds[pos] = 0;
		pos--;
		if (pos==-1)
			break;
		inds[pos]++;
	}
}


/* 
*	Expand the current vector into a new (larger) multi dimensional vector.
 *   inds_in_big gives the indices of the small array variables in the big array
 *
 */
void MulDimArr::ExpandAndAdd(MulDimArr & big_to_add_to, vector<int> & inds_of_small_in_big)
{
	int vi;
	// TODO: Verify that var_sizes_big agrees with the var_sizes of the smaller (this) array
	//	MulDimArr big_arr(var_sizes_big);
	// This will index the big array
	vector<int> inds_for_big;
	
	// Initialize vector of indices for big array to zeros. 
	for (int i=0; i<big_to_add_to.m_base_sizes.size(); i++)
		inds_for_big.push_back(0);
	
	// We now go over the big indices one by one, and for each, take the value from the indices
	// corresponding to the small (this) array
	double *big_dat = &big_to_add_to.m_dat[0];
	for (vi=0; vi<big_to_add_to.m_n_prodsize; vi++) 
	{
		// Get the flat index in the small array
		int ind = GetFlatIndFromBigSpecial(inds_for_big,inds_of_small_in_big);
		// Put value in big array (note that the "safe" way to do this would be to translate inds_for_big into a flat index, but
		// we don't do this to save computation
		//		big_to_add_to.m_dat[vi] += m_dat[ind];
		(*big_dat)+=m_dat[ind];
		big_dat++;
		// Move to the next indices for big (note again that we assume this corresponds to a flat index of vi+1)
		big_to_add_to.BaseIncSpecial(inds_for_big);
	}
}	


/* 
*	Expand the current vector into a new (larger) multi dimensional vector.
 *   inds_in_big gives the indices of the small array variables in the big array
 *
 */
MulDimArr MulDimArr::Expand(vector<int> & var_sizes_big, vector<int> & inds_of_small_in_big)
{
	int vi;
	// TODO: Verify that var_sizes_big agrees with the var_sizes of the smaller (this) array
	MulDimArr big_arr(var_sizes_big);
	// This will index the big array
	vector<int> inds_for_big;
	
	// Initialize vector of indices for big array to zeros. 
	for (int i=0; i<var_sizes_big.size(); i++)
		inds_for_big.push_back(0);
	
	// We now go over the big indices one by one, and for each, take the value from the indices
	// corresponding to the small (this) array
	for (vi=0; vi<big_arr.m_n_prodsize; vi++) 
	{
		// Get the flat index in the small array
		int ind = GetFlatIndFromBigSpecial(inds_for_big,inds_of_small_in_big);
		// Put value in big array (note that the "safe" way to do this would be to translate inds_for_big into a flat index, but
		// we don't do this to save computation
		big_arr.m_dat[vi] = m_dat[ind];
		// Move to the next indices for big (note again that we assume this corresponds to a flat index of vi+1)
		big_arr.BaseIncSpecial(inds_for_big);
	}
	
	return big_arr;
}	


MulDimArr & MulDimArr::operator*=(double val)
{
	double *p;
	for (p=m_dat; p<m_ep; p++)
		(*p)*=val;
	return (*this);
}

MulDimArr & MulDimArr::operator=(double val)
{
	double *p;
	for (p=m_dat; p<m_ep; p++)
		(*p)=val;
	return (*this);
}

MulDimArr & MulDimArr::operator+=(MulDimArr & v)
{
	double *p,*pv;
	for (p=m_dat, pv=v.m_dat; p<m_ep; p++, pv++)
		(*p)+=(*pv);
	return (*this);
}

MulDimArr & MulDimArr::operator-=(MulDimArr & v)
{
	double *p,*pv;
	for (p=m_dat, pv=v.m_dat; p<m_ep; p++, pv++)
		(*p)-=(*pv);
	return (*this);
}

double MulDimArr::Max(int &max_at)
{
	double m = m_dat[0];
	max_at = 0;
	for (int i=1;i<m_n_prodsize; i++) {
		if (m_dat[i]>=m) 
		{
			max_at = i;
			m=m_dat[i];
		}
	}
	return m;
}

// For a given subset of variables, for any assignment to the subset maximize over
// the value of all variables outside the subset
// TODO: Do this for multiple subsets simultaneously. Should be much more efficient!!
void MulDimArr::max_into_multiple_subsets_special(vector<vector<int> > & all_subset_inds, vector<MulDimArr> & all_max_res)
{
	// Prepare a MulDimArr into which we'll maximize
	// First set its size
	int i,si;
	int nSubsets = all_subset_inds.size();

	bool *b_need_max = new bool[nSubsets];
	
	for (si=0; si<nSubsets; si++)
	{
		// If the subset equals the big array	then maximizing would give us the subset (assuming there is no reordering, which I think we can)
		if (all_subset_inds[si].size()==m_base_sizes.size())
		{
			all_max_res[si] = (*this);
			b_need_max[si] = 0;
		}
		else	
		{
			all_max_res[si] = -1e9;
			b_need_max[si] = 1;
		}			
	}
	
	// Go over all values of the big array (this). For each check if its
	// value on the subset is larger than the current max 
	// NOTE: We make the same assumption here as in the Expand function, regaring
	// the correspondence between the flat index (vi) and the the index vector inds_for_big
	vector<int> inds_for_big;
	for (i=0; i< m_base_sizes.size(); i++)
		inds_for_big.push_back(0);

	for (int vi=0; vi<m_n_prodsize; vi++)
	{
		for (si=0; si<nSubsets; si++)
		{
			if (!b_need_max[si])
				continue;
			int flat_subind = all_max_res[si].GetFlatIndFromBigSpecial(inds_for_big,all_subset_inds[si]);
			all_max_res[si][flat_subind] = max(all_max_res[si][flat_subind],m_dat[vi]);
		}
		BaseIncSpecial(inds_for_big);
	}

	delete [] b_need_max;
}

void MulDimArr::Write(ofstream & ofs)
{
	for (int i=0; i<m_n_prodsize; i++)
		ofs << m_dat[i] << " ";
	ofs << endl;
}

void MulDimArr::Read(ifstream & ifs)
{
	for (int i=0; i<m_n_prodsize; i++)
		ifs >> m_dat[i];
}
