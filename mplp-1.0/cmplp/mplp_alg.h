/*
 *  mplp_alg.h
 *  mplp
 *
 *  Created by Amir Globerson and David Sontag on 6/25/08.
 *  Copyright 2008 MIT. All rights reserved.
 *
 */
using namespace std;
#include <vector>
#include <iostream>
#include <map>

#include "muldim_arr.h"

class Region
{
public:
	// The variable size for each variable in the region
//	Vec m_var_sizes;
	// The size of the state-space for the region
	int m_region_prodsize;
	// The variables in the region
	vector<int> m_region_inds;
//	vector<int> m_var_sizes;	
	MulDimArr m_region_lambda;
	// Specifies the indices of the intersection sets
	vector<int> m_intersect_inds;
	// Specifies the indices corresponding to the position of each
	// intersection set in the region
	vector<vector<int> > m_inds_of_intersects;
	// Contains the messages from each region to its intersection sets
	vector<MulDimArr> m_msgs_from_region;
	vector<int> m_var_sizes; 
	
	Region(vector<int> & region_inds, MulDimArr & region_lambda, vector<vector<int> > & all_intersects, vector<int> & intersect_inds, vector<int> & var_sizes);
	double GetValFromGlobal(vector<int> & global_assignment);
	
	void UpdateMsgs(vector<MulDimArr> & sum_into_intersects);
	int Get_nVars() {return m_var_sizes.size();};
};

class MPLPAlg
{
public:
	vector<vector<int> > m_all_intersects;
	vector<Region> m_all_regions;
	vector<vector<int> > m_all_region_inds;	
	vector<MulDimArr> m_sum_into_intersects;
	vector<MulDimArr> m_all_lambdas;
	vector<vector<int> > m_all_region_intersects;
	vector<int> m_var_sizes;			
	vector<int> m_decoded_res;
	vector<double> m_objhist;
	vector<double> m_inthist;
	vector<double> m_timehist;
	
	// This map allows us to quickly look up the index of edge intersection sets
	map<pair<int, int>, int> m_intersect_map;
	
	MPLPAlg(char *regions_fname, char *intersects_fname, char *region_intersects_fname, char *lambdas_fname, char *varsizes_fname);
	void RunMPLP(int niter, double obj_del_thr, double int_gap_thr);
	
	// Add a new region and return its index. intersect_inds refers to the index of the intersection sets that this
	// region intersects with (that is, the index into m_all_intersects)
	int AddRegion(vector<int> & inds_of_vars, MulDimArr & region_lambda,  vector<int> & intersect_inds);
	// As before, but use all zero potential
	int AddRegion(vector<int> & inds_of_vars, vector<int> & intersect_inds);
	// As in the Matlab code, for now we will assume that an intersetion set is added before adding the regions that
	// intersect with it
	int AddIntersectionSet(vector<int> & inds_of_vars);
	// Find the index number into m_all_intersects of a given set of variables' intersection set.
	// Returns -1 if not found.
	int FindIntersectionSet(vector<int> & inds_of_vars);
	
	void Read(char *fname);
	void Write(char *msgs_fname,char *res_fname, char *suminto_fname, char *objhist_fname, char *inthist_fname, char *timehist_fname);
};
