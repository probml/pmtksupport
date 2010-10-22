/*
 *  mplp_alg.cpp
 *  mplp
 *
 *  Created by Amir Globerson and David Sontag on 6/25/08.
 *  Copyright 2008 MIT. All rights reserved.
 *
 */

#include "mplp_alg.h"
#include <iostream>
#include <sstream>
#include <fstream>
#include <ctime>
#include <algorithm>

void print_int_vec(vector<int> v)
{
	for (int i=0; i<v.size(); i++)
		cout << v[i] << " ";
//	cout << endl;
}

void print_double_vec(vector<double> v)
{
	for (int i=0; i<v.size(); i++)
		cout << v[i] << " ";
	cout << endl;
}

Region::Region(vector<int> & region_inds, MulDimArr & region_lambda, vector<vector<int> > & all_intersects, vector<int> & intersect_inds, vector<int> & var_sizes): m_region_lambda(region_lambda), m_region_inds(region_inds), m_intersect_inds(intersect_inds)
{
	// Find the indices of each intersection within the region. Also intialize the message into that intersection
	for (int si=0; si<m_intersect_inds.size(); si++)
    {
		vector<int> tmp_inds_of_intersects;
		
		vector<int> curr_intersect = all_intersects[m_intersect_inds[si]];
		
		vector<int> intersect_var_sizes;
		// Go over all variables in the intersection set
		for (int i=0; i< curr_intersect.size(); i++)
		{
			int var_in_intersect = curr_intersect[i];
			intersect_var_sizes.push_back(var_sizes[var_in_intersect]);
			
			vector<int>::iterator iter = find(m_region_inds.begin(), m_region_inds.end(), var_in_intersect);
			

			// Verify that the current intersection variable indeed appears in this region, and get the index where it appears
			if (iter == m_region_inds.end())
			{
				cerr << "Intersection set contains variables not in region" << endl;
				return;
			}
			else
				tmp_inds_of_intersects.push_back(iter-m_region_inds.begin());
		}
		m_inds_of_intersects.push_back(tmp_inds_of_intersects);
		
		// This will initialize the message and  set it to zero
		MulDimArr curr_msg(intersect_var_sizes);
		curr_msg = 0;
		
		m_msgs_from_region.push_back(curr_msg);
    }
	// Calculate the size of the region state space (although this should already be in the lambda object, so we should
    // probably avoid this multiplicity)
    for (int i=0; i<region_inds.size(); i++) {
		int curr_size = var_sizes[region_inds[i]];
        m_var_sizes.push_back(curr_size);
	}
}

double Region::GetValFromGlobal(vector<int> & global_assignment)
{
	// global_assignment gives values for all variables. We need to take
	// the values for the variables corresponding to the current region,
	// put them into a vector<int>, and get the value of the region's
	// lambda at the relevant indices
	// Check if lambda is an empty vector, in which case return 0
	if (m_region_lambda.m_n_prodsize==0)
	  return 0;
	// Prepare vector with indices in Region
	vector<int> indices;
	
	for (int vi=0; vi<m_region_inds.size(); vi++)
		indices.push_back(global_assignment[m_region_inds[vi]]);
	return m_region_lambda.GetVal(indices);
}

void Region::UpdateMsgs(vector<MulDimArr> & sum_into_intersects)
{

	/* First do the expansion:
	1. Take out the message into the intersection set from the current cluster
	2. Expand it to the size of the region
	3. Add this for all intersection sets
	*/
	MulDimArr expansion(m_var_sizes);
	
	// The following if is for the case where lambda is empty (as we use for the cluster potentials) which means we
	// should set it to zero.
	if (m_region_lambda.m_n_prodsize==0)
		expansion = 0;
	else
		expansion = m_region_lambda;

	// Will store the total messages going into the intersection, but not from the Region
	vector<MulDimArr> lam_minus_region;
	for (int si=0; si<m_intersect_inds.size(); si++)
    {
		int curr_intersect = m_intersect_inds[si];
		MulDimArr tmp(sum_into_intersects[curr_intersect]);
		tmp-=m_msgs_from_region[si];
		lam_minus_region.push_back(tmp);

		vector<int> & curr_inds_of_intersect = m_inds_of_intersects[si];
		// If the intersection has the same size as the region, we assume they are the same, and therefore
		// no need to expand 
		// NOTE NOTE: There will be problems with the above if the intersection has the same indices
		// but rearranged.
		if (Get_nVars()==curr_inds_of_intersect.size())
			expansion+= tmp;
		else
		{
			tmp.ExpandAndAdd(expansion,curr_inds_of_intersect);
			//			expansion+= tmp_expanded;  
		}
		// Take out the current incoming message
		sum_into_intersects[curr_intersect]-= m_msgs_from_region[si];	
    }
	// Update messages
	expansion.max_into_multiple_subsets_special(m_inds_of_intersects,m_msgs_from_region);
	int sC = m_intersect_inds.size();
	for (int si=0; si<m_intersect_inds.size(); si++)
    {
		// Take out previous message
		int curr_intersect = m_intersect_inds[si];		
		MulDimArr curr_max_into_subset;
		// Update message
		m_msgs_from_region[si]*= 1.0/sC;
		m_msgs_from_region[si] -= lam_minus_region[si];	

		// Put in current message
		sum_into_intersects[curr_intersect]+= m_msgs_from_region[si];
    }	
	return;
}

MPLPAlg::MPLPAlg(char *regions_fname, char *intersects_fname, char *region_intersects_fname, char *lambdas_fname, char *varsizes_fname)
{
	ifstream regions_infile(regions_fname, ios_base::in);
	ifstream intersects_infile(intersects_fname, ios_base::in);
	ifstream regions_intersects_infile(region_intersects_fname, ios_base::in);  
	ifstream lambdas_infile(lambdas_fname, ios_base::in);	
	ifstream varsizes_infile(varsizes_fname, ios_base::in);		
	
	string line;
	vector<vector<double> > all_lambdas;
	
	/* IMPORTANT NOTE: We assume here the input indices are one based (so we subtract one to translate to C indices */
	/* TODO: We keep the initialization information in memebers of MPLPAlg. This is not necessary, since the information
	   is also kept in the regions themselves, where it is actually used. Need to remove redundant members */
	   
	// Read the regions from the region file (the i'th row gives the indices of the variables in the i'th region 
	while (! regions_infile.eof() )
    {
		getline(regions_infile, line);
		istringstream str(line);
		int foo;
		vector<int> curr_region;
		while (str>>foo)
			curr_region.push_back(foo-1);
		if (curr_region.size()==0)
			break;
		m_all_region_inds.push_back(curr_region);
    }	

	cout << "Read " << m_all_region_inds.size() << " regions." << endl;
	// Read the intersection sets (the i'th row gives the indices of variables in the i'th intersection sets 	
	while (! intersects_infile.eof() )
    {
		getline(intersects_infile, line);
		istringstream str(line);
		int foo;
		vector<int> curr_intersect;
		while (str>>foo)
			curr_intersect.push_back(foo-1);
		if (curr_intersect.size()==0)
			break;
				
		m_all_intersects.push_back(curr_intersect);
		
		// If this is an edge, insert into the map
		if(curr_intersect.size() == 2)
		{
			// First sort
			vector<int> tmp_inds(curr_intersect);
			sort(tmp_inds.begin(), tmp_inds.end());
			
			// Then insert
			m_intersect_map.insert(pair<pair<int,int>,int>(pair<int,int>(tmp_inds[0], tmp_inds[1]), m_all_intersects.size()-1));
		}
		
    }	

	cout << "Read " << m_all_intersects.size() << " intersection sets." << endl;

	// Read which intersects belong to which region (the i'th row gives the indices of intersection sets that belong to the i'th region 		
	while (! regions_intersects_infile.eof() )
    {
		getline(regions_intersects_infile, line);
		istringstream str(line);
		int foo;
		vector<int> curr_intersect;
		while (str>>foo)
			curr_intersect.push_back(foo-1);
		if (curr_intersect.size()==0)
			break;
		m_all_region_intersects.push_back(curr_intersect);
    }	

	cout << "Read " << m_all_region_intersects.size() << " region intersection inputs." << endl;
	
	// Read the sizes of the variables. This is just a vector
	while (! varsizes_infile.eof() )
    {	
		int foo;
		
		if (varsizes_infile >> foo)
			m_var_sizes.push_back(foo);
    }
	cout << "There are " << m_var_sizes.size() << " variables in the model." << endl;
	
	
	// Read the potentials. Put each one in a vector<double> and then transform into MulDimArr
	while (! lambdas_infile.eof() )	
    {
		getline(lambdas_infile, line);
		istringstream str(line);
		double foo;
		vector<double> curr_lambda;
		while (str>>foo)
			curr_lambda.push_back(foo);
		// If it's an empty line, it could be just a lambda that is zeros (unless we have read all the regions)
		if (curr_lambda.size()==0 && (all_lambdas.size()==m_all_region_inds.size()))
			break;
		all_lambdas.push_back(curr_lambda);
    }

	cout << "Read " << all_lambdas.size() << " lambda vectors." << endl;
	
	// Initialize sum into intersectios. 
	// NOTE: As written, this assumes all messages are initialized to zero 
	for (int si=0; si<m_all_intersects.size(); si++)
    {
		// Figure out the variable sizes for the intersection (put in subset_size)
		vector<int> curr_intersection_set = m_all_intersects[si];
		vector<int> subset_size ;
		for (int i=0; i<curr_intersection_set.size(); i++)
			subset_size.push_back(m_var_sizes[curr_intersection_set[i]]);
		MulDimArr curr_suminto(subset_size);
		curr_suminto = 0;
		m_sum_into_intersects.push_back(curr_suminto);
    }

	cout << "Initializing regions" << endl;
	clock_t start,finish;
	double time;

	start = clock();

	// Initialize regions
	for (int ri=0; ri<m_all_region_inds.size(); ri++)
	{
	  //	       cout << "Init region " << ri << endl; 
		// Create size for this region
		vector<int> region_var_sizes;
		int i;

		  for (i=0; i< m_all_region_inds[ri].size(); i++)
		    region_var_sizes.push_back(m_var_sizes[m_all_region_inds[ri][i]]);
		  
		MulDimArr *curr_lambda;

		if (all_lambdas[ri].size()!=0)
		{
		  curr_lambda = new MulDimArr(region_var_sizes);
		  // Assume all_lambdas is given as a "flat" vector and put it into curr_lambda
		  for (i=0; i< curr_lambda->m_n_prodsize; i++)
		    (*curr_lambda)[i] = all_lambdas[ri][i];
		}
		else
		  curr_lambda = new MulDimArr; // Empty constructor

		// Initialize current region
		Region curr_region(m_all_region_inds[ri],*curr_lambda,m_all_intersects, m_all_region_intersects[ri], m_var_sizes);
		delete curr_lambda;
		m_all_regions.push_back(curr_region);
	}
	
	finish = clock();
	time = (double(finish)-double(start))/CLOCKS_PER_SEC;	
	cout << "Region initialization " <<  " took " << time << endl;

	// Initialize output vector
	for (int i=0; i<m_var_sizes.size(); i++)
		m_decoded_res.push_back(0);
}

void MPLPAlg::RunMPLP(int niter, double obj_del_thr, double int_gap_thr)
{
	int ri;
	double last_obj=1e10;
	
	cout << "Beginning message passing" << endl;
	

	// Next couple of lines will pass the message from the single node to the single node intersection set (basically
	// just pass the evidence). This is not strictly necessary and I only do it to be compatible with the way it's done in the matlab code

	for (ri=0; ri<m_all_regions.size(); ri++)
	{
		if (m_all_regions[ri].Get_nVars()==1)
			m_all_regions[ri].UpdateMsgs(m_sum_into_intersects);
	}


	// Perform the GMPLP updates as in the NIPS07 paper
	for (int it=0; it<niter; it++)
    {
		clock_t start,finish;
		double time;

		start = clock();
		for (ri=0; ri<m_all_regions.size(); ri++)
		{
			// For the single node regions, there is a need to pass a message
			// only at the first iteration (since it will pass the local potential,
			// and will be the same for all iterations).
		  if ((m_all_regions[ri].Get_nVars()==1) && (it>0))
		  	continue;
		  m_all_regions[ri].UpdateMsgs(m_sum_into_intersects);
		}

		// Calculate objective
		double obj=0;
		for (int si=0; si<m_sum_into_intersects.size(); si++)
		{
			int max_at;
			obj+= m_sum_into_intersects[si].Max(max_at);
			// If this is a singleton, keep its value (so that we also have an integral assignment). 
			// NOTE: Here we assume that all singletons are intersection sets. Otherwise, some variables will not be decoded here
			if (m_all_intersects[si].size()==1)
				m_decoded_res[m_all_intersects[si][0]] = max_at;
		}
		// Get value of decoded solution
		double int_val=0;
		for (ri=0; ri<m_all_regions.size(); ri++)
			int_val+=m_all_regions[ri].GetValFromGlobal(m_decoded_res);
		// Print results
		finish = clock();
		time = (double(finish)-double(start))/CLOCKS_PER_SEC;	

		double int_gap = obj-int_val;
		double obj_del = last_obj-obj;
		m_objhist.push_back(obj);
		m_inthist.push_back(int_val);
		m_timehist.push_back(time);

		cout << "Iter=" << (it+1) << " Time=" << time << " Objective=" << obj <<  " Decoded=" << int_val << " ObjDel=" <<  obj_del << " IntGap=" << int_gap << endl;
		if (obj_del<obj_del_thr)
			break;
		if (int_gap<int_gap_thr)
			break; 
		last_obj = obj;
    }
	return;
}

void MPLPAlg::Write(char *msgs_fname,char *res_fname, char *suminto_fname, char *objhist_fname, char *inthist_fname, char *timehist_fname)
{
	// Write the messages themselves
	ofstream ofs_messages(msgs_fname);
	ofstream ofs_suminto(suminto_fname);
	ofstream ofs_objhist(objhist_fname);
	ofstream ofs_inthist(inthist_fname);
	ofstream ofs_timehist(timehist_fname);
	ofstream ofs_res(res_fname);
	
	int ri,vi,si,it;
	
	for (ri=0; ri<m_all_regions.size(); ri++)
		for (si=0; si< m_all_regions[ri].m_msgs_from_region.size(); si++)
			m_all_regions[ri].m_msgs_from_region[si].Write(ofs_messages);
	ofs_messages.close();
	
	
	for (vi=0; vi< m_var_sizes.size(); vi++)
		ofs_res << m_decoded_res[vi] << " ";
	ofs_res << endl;
	ofs_res.close();
	
	for (si=0; si<m_sum_into_intersects.size(); si++)
		m_sum_into_intersects[si].Write(ofs_suminto);
	
	for (it=0; it<m_objhist.size(); it++)
		ofs_objhist << m_objhist[it] << " ";
	ofs_objhist << endl;
	
	for (it=0; it<m_inthist.size(); it++)
		ofs_inthist << m_inthist[it] << " ";
	ofs_inthist << endl;

	for (it=0; it<m_timehist.size(); it++)
		ofs_timehist << m_timehist[it] << " ";
	ofs_timehist << endl;
	
	ofs_messages.close();
	ofs_suminto.close();
	ofs_objhist.close();
	ofs_inthist.close();
	ofs_res.close();
	ofs_timehist.close();
}

void MPLPAlg::Read(char *fname)
{
	// Read messages from file, and calculate sum into intersections. Assume sum_into_intersections has been set
	// to zero (which it is in the constructor)
	ifstream ifs_messages(fname);
	if (!ifs_messages.is_open())
	{
		cout << "No messages file" << endl;
		return;
	}	
	for (int ri=0; ri<m_all_regions.size(); ri++)
		for (int si=0; si< m_all_regions[ri].m_msgs_from_region.size(); si++)
		{
			m_all_regions[ri].m_msgs_from_region[si].Read(ifs_messages);
			m_sum_into_intersects[m_all_regions[ri].m_intersect_inds[si]]+=m_all_regions[ri].m_msgs_from_region[si];
		}	
	ifs_messages.close();
}


int MPLPAlg::AddRegion(vector<int> & inds_of_vars, MulDimArr & region_lambda,  vector<int> & intersect_inds)
{
	// This will also initialize the messages to zero, which is what we want
	Region new_region(inds_of_vars, region_lambda, m_all_intersects, intersect_inds, m_var_sizes);
	m_all_regions.push_back(new_region);
	return m_all_regions.size()-1;
}

int MPLPAlg::AddRegion(vector<int> & inds_of_vars, vector<int> & intersect_inds)
{
    // initialize parameters to all zero.
	
	// Calculate the sizes of the variables in this set
	vector<int> sizes;
	for (int i=0; i< inds_of_vars.size(); i++)
		sizes.push_back(m_var_sizes[inds_of_vars[i]]);
	
	MulDimArr region_lambda(sizes);
	return AddRegion(inds_of_vars, region_lambda, intersect_inds);
}

int MPLPAlg::AddIntersectionSet(vector<int> & inds_of_vars)
{	
	m_all_intersects.push_back(inds_of_vars);
	// Calculate the sizes of the variables in this set
	vector<int> sizes;
	for (int i=0; i< inds_of_vars.size(); i++)
		sizes.push_back(m_var_sizes[inds_of_vars[i]]);
	
	// If this is an edge, insert into the map
	if(inds_of_vars.size() == 2)
	{
		// First sort
		vector<int> tmp_inds(inds_of_vars);
		sort(tmp_inds.begin(), tmp_inds.end());

		// Then insert
		m_intersect_map.insert(pair<pair<int,int>,int>(pair<int,int>(tmp_inds[0], tmp_inds[1]), m_all_intersects.size()-1));
	}
	
	MulDimArr new_arr(sizes);
	new_arr = 0;
	m_sum_into_intersects.push_back(new_arr);
	return m_all_intersects.size()-1;
}

int MPLPAlg::FindIntersectionSet(vector<int> & inds_of_vars)
{
	// Sort the indices to make lookup and comparison easy
	vector<int> tmp_inds_of_vars(inds_of_vars);
	sort(tmp_inds_of_vars.begin(), tmp_inds_of_vars.end());

	// Is this an edge? If so, do map lookup
	if(tmp_inds_of_vars.size() == 2)
	{
		map<pair<int,int>, int>::iterator iter = m_intersect_map.find(pair<int,int>(tmp_inds_of_vars[0], tmp_inds_of_vars[1]));
        if( iter != m_intersect_map.end() ) // If the edge is found
			return iter->second;
        else {
			cerr << "FindIntersectionSet asked to find an _edge_ intersection set that does not exist!" << endl;
			return -1;
		}
	}
	
	for(int i=0; i < m_all_intersects.size(); i++)
	{
		// copy, then sort
		vector<int> tmp_inds(m_all_intersects[i]);
		sort(tmp_inds.begin(), tmp_inds.end());
		
		if(tmp_inds == tmp_inds_of_vars)
			return i;
	}
	
	cerr << "FindIntersectionSet asked to find an intersection set that does not exist!" << endl;
	return -1;
}
