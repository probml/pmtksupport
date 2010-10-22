/*
 *  muldim_arr.h
 *  mplp
 *
 *  Created by Amir Globerson on 6/24/08.
 *  Copyright 2008 MIT. All rights reserved.
 *
 */
using namespace std;
#include <vector>
#include <iostream>
#include "matrix.h"
#include <string.h> 
#ifndef _MULDIMARR
#define _MULDIMARR

#define BASE2DEC_PAIR(x1,x2,base1,base2) (x1*base2+x2)

class MulDimArr 
{
public:
	vector<int> m_base_sizes;
	int m_n_prodsize;
	double *m_dat;
	double *m_ep;
	
	// Initialize to all zero
	MulDimArr(vector<int> & base_sizes); 
	
	// Copy constructor
	MulDimArr(const MulDimArr & v);
	
	MulDimArr() 
	{
		m_n_prodsize = 0;
		m_dat = NULL;
	};
	~MulDimArr()
	{
		if (m_dat!=NULL)
			delete [] m_dat;
	}
	
	
	MulDimArr & operator=(const MulDimArr & v);
	MulDimArr & operator=(double val);
	MulDimArr & operator*=(double val);
	MulDimArr & operator+=(MulDimArr & v);
	MulDimArr & operator-=(MulDimArr & v);
	double & operator[](int i) {return m_dat[i];}
	double Max(int &max_at);
	void print();
	void print_with_inds();
	void Write(ofstream & ofs);
	void Read(ifstream & ifs);
	
	inline int GetFlatInd(vector<int> base_inds);
	inline int GetFlatIndFromBig(vector<int> big_base_inds, vector<int> inds_in_big);
	int GetFlatIndFromBigSpecial(vector<int> & big_base_inds, vector<int> & inds_in_big);

	double GetVal(vector<int> & indices);
	MulDimArr Expand(vector<int> & var_sizes_big, vector<int> & inds_in_big);
	void ExpandAndAdd(MulDimArr & big_to_add_to, vector<int> & inds_of_small_in_big);

	inline void BaseInc(vector<int> & inds);
	inline void BaseIncSpecial(vector<int> & inds);

	void  max_into_multiple_subsets_special(vector<vector<int> > & all_subset_inds, vector<MulDimArr> & all_maxes);
};

// Return the index in the flat multi-dimensional array corresponding to
// the given multi-index
// NOTE: Since this function is called to get an intersection index out of a big index,
// as long as we don't use intersections of more than two, we only need the cases given here.
inline int MulDimArr::GetFlatIndFromBigSpecial(vector<int> & big_base_inds, vector<int> & inds_in_big)
{
	int y,ind1,ind2;
	
	switch (inds_in_big.size())
	{
		case 1:
			y = big_base_inds[inds_in_big[0]];
			break;
		case 2:
//			cout << "Not good here too" << endl;
			ind1 = inds_in_big[0];
			ind2 = inds_in_big[1];
			y = BASE2DEC_PAIR(big_base_inds[ind1],big_base_inds[ind2],m_base_sizes[0],m_base_sizes[1]);
			break;
		default:
			cout << "GetFlatIndFromBigSpecial problem" << endl;
			break;
	}
	return y;
}

inline void MulDimArr::BaseIncSpecial(vector<int> & inds)
{
	switch (inds.size())
	{
		case 1:
			inds[0]++;
			break;
		case 2:
			inds[1]++;
			if (inds[1]==m_base_sizes[1])
			{
				inds[1] =0;
                inds[0]++;
            }	
			break;	
		case 3:
			inds[2]++;
			if (inds[2]==m_base_sizes[2])
			{
				inds[2]=0;
				inds[1]++;
				if (inds[1]==m_base_sizes[1])
				{
					inds[1]=0;
					inds[0]++;
				}
			}
			break;
		default: 
			cout << "Not supported yet" << endl;
	}
}

#endif
