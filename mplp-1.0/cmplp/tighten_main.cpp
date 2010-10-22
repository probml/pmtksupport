/*
 *  tighten_main.cpp
 *  mplp
 *
 *  Created by Amir Globerson and David Sontag on 8/10/08.
 *  Copyright 2008 MIT. All rights reserved.
 *
 */
#include <iostream>

using namespace std;
#include <vector>
#include <ctime>
#include <algorithm>
#include <set>
#include "mplp_alg.h"
#define MAX_TIGHT_ITERS 1000000
#define GRAPH_TYPE_DEFAULT 0
#define GRAPH_TYPE_GRID 1
#define sub2ind(numCols, numRows, c, r) ((r)*numCols + c)
#define Inf 9999999999.9
typedef map<pair<int,int>, int> mapType;


struct TripletCluster {
	double bound;
	int i,j,k;
	int ij_intersect_loc, jk_intersect_loc, ki_intersect_loc;
	
	bool operator <(const TripletCluster & rhs) const {
		return bound < rhs.bound;
	}
};

struct SquareCluster {
	double bound;
	int i,j,k,l;
	int ij_intersect_loc, kj_intersect_loc, kl_intersect_loc, il_intersect_loc;
	
	bool operator <(const SquareCluster & rhs) const {
		return bound < rhs.bound;
	}
};

double maximizeIndependently(vector<MulDimArr*> & beliefs)
{
	double sum=0.0;
	int max_at; // not actually needed
	for(int i=0; i < beliefs.size(); i++)
	{
		sum += beliefs[i]->Max(max_at);
	}
	
	return sum;
}

double getValCycle(vector<MulDimArr*> & beliefs, vector<bool> & b_transpose, vector<int> & assignments)
{
	double sum=0.0;
	vector<int> inds; inds.push_back(-1); inds.push_back(-1); // temp
	
	// All except the last edge
	for(int i=0; i < beliefs.size()-1; i++)
	{
		inds[b_transpose[i]?1:0] = assignments[i];
		inds[b_transpose[i]?0:1] = assignments[i+1];
		sum += beliefs[i]->GetVal(inds);
	}
	
	// Now do last edge
	inds[b_transpose[beliefs.size()-1]?1:0] = assignments[beliefs.size()-1];
	inds[b_transpose[beliefs.size()-1]?0:1] = assignments[0];
	sum += beliefs[beliefs.size()-1]->GetVal(inds);
	
	return sum;
}

double maximizeCycle(vector<MulDimArr*> & beliefs, vector<bool> & b_transpose)
{
	double max_val = -Inf;
	
	// Fix value of the first variable
	int first_var_size = beliefs[0]->m_base_sizes[b_transpose[0]?1:0];
	int second_var_size = beliefs[0]->m_base_sizes[b_transpose[0]?0:1];
	for(int vo=0; vo < first_var_size; vo++)
	{
		vector<int> inds; inds.push_back(-1); inds.push_back(-1); // temp
		inds[b_transpose[0]?1:0] = vo;
		
		// Do first edge (construct initial field)
		vector<double> field;
		for(int v2=0; v2 < second_var_size; v2++)
		{
			inds[b_transpose[0]?0:1] = v2;
			field.push_back(beliefs[0]->GetVal(inds));
		}
		
		// Go over rest of edges, except last (which has to be treated specially)
		for(int i=1; i < beliefs.size()-1; i++)
		{
			vector<double> new_field;
			for(int v2=0; v2 < beliefs[i]->m_base_sizes[b_transpose[i]?0:1]; v2++)
			{
				inds.clear(); inds.push_back(-1); inds.push_back(-1); // temp
				inds[b_transpose[i]?0:1] = v2;
				
				// Take max
				double tmp_max_val = -Inf;
				for(int v1=0; v1 < field.size(); v1++)
				{
					inds[b_transpose[i]?1:0] = v1;
					tmp_max_val = max(tmp_max_val, field[v1]+beliefs[i]->GetVal(inds));
				}
				new_field.push_back(tmp_max_val);
			}
			field.clear(); // necessary?
			field = new_field;
		}
		
		// Do last edge (fix endpoint value to vo)
		inds.clear(); inds.push_back(-1); inds.push_back(-1); // temp
		inds[b_transpose[b_transpose.size()-1]?0:1] = vo;

		// Take max
		double tmp_max_val = -Inf;		
		for(int v1=0; v1 < field.size(); v1++)
		{
			inds[b_transpose[b_transpose.size()-1]?1:0] = v1;
			tmp_max_val = max(tmp_max_val, field[v1]+beliefs[beliefs.size()-1]->GetVal(inds));
		}
		
		max_val = max(max_val, tmp_max_val);
	}

	return max_val;
}


int main( int argc, char *argv[] )
{
	int niter = 10;
	int niter_later = 10;
	int nclus_to_add = 1;
	int graph_type = 0;
	int numRows, numCols;
	double obj_del_thr;
	double int_gap_thr;
	clock_t start,finish;
	double time;
	
	if (argc<7)
	{
		printf("Syntax: algo_triplet niter niter_later nclus_to_add obj_del_thr int_gap_thr graph_type (num_rows num_cols)\n");
		printf("        (graph_type: default %d, grid %d)\n", GRAPH_TYPE_DEFAULT, GRAPH_TYPE_GRID);
		return 0;
	}
	
	sscanf(argv[1],"%d",&niter);
	sscanf(argv[2],"%d",&niter_later);
	sscanf(argv[3],"%d",&nclus_to_add);
	sscanf(argv[4],"%lg",&obj_del_thr);
	sscanf(argv[5],"%lg",&int_gap_thr);	
	sscanf(argv[6],"%d",&graph_type);
	
	if(graph_type == GRAPH_TYPE_GRID)
	{
		// Needs to have additional arguments for num_rows and num_cols
		if (argc<9)
		{
			printf("*** Since you say this is a grid, can you please provide num_rows and num_cols? ***\n");
			printf("Syntax: cmplp niter niter_later nclus_to_add obj_del_thr int_gap_thr graph_type (num_rows num_cols)\n");
			printf("        (graph_type: default %d, grid %d)\n", GRAPH_TYPE_DEFAULT, GRAPH_TYPE_GRID);
			return 0;
		}

		sscanf(argv[7],"%d",&numRows);
		sscanf(argv[8],"%d",&numCols);
	}
	
	printf("niter=%d\nniter_later=%d\nnclus_to_add=%d\nobj_del_thr=%lg\nint_gap_thr=%lg\ngraph_type=%d\n",niter, niter_later, nclus_to_add,obj_del_thr,int_gap_thr, graph_type);
	
	// Load in the MRF and initialize GMPLP state
	MPLPAlg mplp("regions.txt","intersects.txt","region_intersects.txt","lambdas.txt","var_sizes.txt");

	// Initialize adjacency list (filled in later) TODO: only do this when needed
	vector<int> adjacency_list[mplp.m_var_sizes.size()];
	int nNewClusters;
	
	// TODO: if we want to re-start computation after adding clusters, must also read/write regions etc.
	//	mplp.Read("msgs.txt");
	
	printf("Initially running MPLP for %d iterations\n", niter);
	mplp.RunMPLP(niter,obj_del_thr,int_gap_thr);

	for(int i=1; i<MAX_TIGHT_ITERS; i++)  // Break when problem is solved
	{
		printf("\n\nOuter loop iteration %d\n----------------------\n", i);

		// printf("Writing intermediate results to output files...\n");
		// mplp.Write("msgs.txt","res.txt","suminto.txt","objhist.txt","inthist.txt","timehist.txt");

		// Is problem solved? If so, break.
		double best_decoding = *max_element(mplp.m_inthist.begin(), mplp.m_inthist.end());  // TODO opt: keep track of max
		double int_gap = mplp.m_objhist.back() - best_decoding;
		if(int_gap < int_gap_thr)
		{
			printf("Done! Integrality gap less than %lg\n", int_gap_thr);
			break;
		}

		if(i==1 && graph_type == GRAPH_TYPE_DEFAULT)
		{
			cout << "Doing pre-processing for adding triplet clusters." << endl;
			
			// How many triangles are there in the graph?
			nNewClusters = 0;
			
			// Construct adjacency list for the graph
			// Iterate over all of the edges (we do this by looking at the edge intersection sets)
			for(mapType::const_iterator it = mplp.m_intersect_map.begin(); it != mplp.m_intersect_map.end(); ++it)
			{
				
				// Get the two nodes i & j
				int i=it->first.first; int j=it->first.second;
				adjacency_list[i].push_back(j);
				adjacency_list[j].push_back(i);
			}
			// Sort the adjacency list, for fast intersections later
			for(int i=0; i < sizeof(adjacency_list)/sizeof(vector<int>); i++) 
			{
				sort(adjacency_list[i].begin(), adjacency_list[i].end());
			}
			
			// Count the number of triangles
			vector<int>::iterator intersects_iter_end;
			vector<int> commonNodes(mplp.m_var_sizes.size());
			for(mapType::const_iterator it = mplp.m_intersect_map.begin(); it != mplp.m_intersect_map.end(); ++it)
			{
				
				// Get the two nodes i & j
				int i=it->first.first; int j=it->first.second;

				// Now find all neighbors of both i and j to see where the triangles are
				intersects_iter_end = set_intersection(adjacency_list[i].begin(), adjacency_list[i].end(), adjacency_list[j].begin(), adjacency_list[j].end(), commonNodes.begin());

				for(vector<int>::const_iterator n=commonNodes.begin(); n != intersects_iter_end; ++n)
				{
					// Since a triplet shows up three times as an edge plus
					// a node, we only consider it for the case when n<i and n<j
					if(*n < i && *n < j)
						nNewClusters++;
				}
			}
		}
		
		// Heuristic: when the integrality gap is sufficiently small, allow the algorithm
		// more time to run till convergence
		if(int_gap < 1)
		{
			niter_later = max(niter_later, 600);  // TODO opt: don't hard code
			obj_del_thr = min(obj_del_thr, 1e-5);
			printf("Int gap small, so setting niter_later to %d and obj_del_thr to %lg\n", niter_later, obj_del_thr);
		}
		
		// Tighten LP
		cout << "Now attempting to tighten LP relaxation..." << endl;

		start = clock();
		int nClustersAdded = 0;

		
		if(graph_type == GRAPH_TYPE_DEFAULT)
		{
			cout << "Looking for triangle clusters." << endl;

			// TODO: put this elsewhere so that the space isn't re-allocated continuously?
			// Enumerate over all of the edges
			TripletCluster newCluster[nNewClusters];
			
			int index=0;
			
			// Iterate over all of the edge intersection sets
			vector<int>::iterator intersects_iter_end;
			vector<int> commonNodes(mplp.m_var_sizes.size());
			vector<int> tripAssignment; tripAssignment.push_back(-1); tripAssignment.push_back(-1); tripAssignment.push_back(-1);
			for(mapType::const_iterator it = mplp.m_intersect_map.begin(); it != mplp.m_intersect_map.end(); ++it)
			{
				
				// Get the two nodes i & j, and the edge intersection index
				int i=it->first.first; int j=it->first.second;
				int ij_intersect_loc = it->second;
				
				// Now find all neighbors of both i and j to see where the triangles are
				intersects_iter_end = set_intersection(adjacency_list[i].begin(), adjacency_list[i].end(), adjacency_list[j].begin(), adjacency_list[j].end(), commonNodes.begin());
				
				for(vector<int>::const_iterator n=commonNodes.begin(); n != intersects_iter_end; ++n)
				{
					int k = *n;
					
					// Since a triplet shows up three times as an edge plus
					// a node, we only consider it for the case when k<i and k<j
					if(!(k < i && k < j))
						continue;
					
					newCluster[index].i = i;
					newCluster[index].j = j;
					newCluster[index].k = k;

					// Find the intersection sets for this triangle
					newCluster[index].ij_intersect_loc = ij_intersect_loc;
					vector<int> jk_edge; jk_edge.push_back(newCluster[index].j); jk_edge.push_back(newCluster[index].k);
					newCluster[index].jk_intersect_loc = mplp.FindIntersectionSet(jk_edge);
					vector<int> ki_edge; ki_edge.push_back(newCluster[index].k); ki_edge.push_back(newCluster[index].i);
					newCluster[index].ki_intersect_loc = mplp.FindIntersectionSet(ki_edge);					

					// Construct the beliefs for each edge, which will be maximized below
					vector<MulDimArr*> beliefs;
					vector<bool> b_transpose;
					beliefs.push_back(&mplp.m_sum_into_intersects[newCluster[index].ij_intersect_loc]);
					b_transpose.push_back(mplp.m_all_intersects[newCluster[index].ij_intersect_loc][0] != newCluster[index].i); // i first
					beliefs.push_back(&mplp.m_sum_into_intersects[newCluster[index].jk_intersect_loc]);
					b_transpose.push_back(mplp.m_all_intersects[newCluster[index].jk_intersect_loc][0] != newCluster[index].j); // then 'j'
					beliefs.push_back(&mplp.m_sum_into_intersects[newCluster[index].ki_intersect_loc]);
					b_transpose.push_back(mplp.m_all_intersects[newCluster[index].ki_intersect_loc][0] != newCluster[index].k); // then 'k'

					double bound_indep = maximizeIndependently(beliefs);
					
					// Before doing expensive joint maximization, see if we can quickly find optimal assignment
					tripAssignment[0] = mplp.m_decoded_res[i]; tripAssignment[1] = mplp.m_decoded_res[j]; tripAssignment[2] = mplp.m_decoded_res[k];
					double bound_quick = getValCycle(beliefs, b_transpose, tripAssignment);
					if(bound_indep == bound_quick)
					{
						newCluster[index].bound = 0;
					} else {
						// Do expensive joint maximization
						newCluster[index].bound = bound_indep - maximizeCycle(beliefs, b_transpose);
					}
					
					index++;
				}
			}
						
			// TODO opt: have a class for a cluster, so we can have different types and sort by bound,
			//       choosing best one by bound.
			//       Make the sorting and adding independent of the type of graph...
			
			// Sort the clusters by the bound
			sort(newCluster, newCluster+nNewClusters);
			
			printf(" -- Considered %d clusters, smallest bound %g, largest bound %g\n", nNewClusters, newCluster[nNewClusters-nclus_to_add].bound, newCluster[nNewClusters-1].bound);
			
			// Add the top nclus_to_add clusters to the relaxation
			for(int clusterId = nNewClusters-1; clusterId >= 0 && nClustersAdded < nclus_to_add; clusterId--)
			{
				// TODO: check that these clusters and intersection sets haven't already been added
				
				// Now add cluster ijk
				vector<int> ijk_inds;
				ijk_inds.push_back(newCluster[clusterId].i); ijk_inds.push_back(newCluster[clusterId].j); ijk_inds.push_back(newCluster[clusterId].k);
				
				vector<int> ijk_intersect_inds;
				ijk_intersect_inds.push_back(newCluster[clusterId].ij_intersect_loc);
				ijk_intersect_inds.push_back(newCluster[clusterId].jk_intersect_loc);
				ijk_intersect_inds.push_back(newCluster[clusterId].ki_intersect_loc);
				
				mplp.AddRegion(ijk_inds, ijk_intersect_inds);
				
				// TODO: log which clusters are chosen...
				
				nClustersAdded++;
			}

		}
		
		
		if(graph_type == GRAPH_TYPE_GRID)
		{
			cout << "Looking for square clusters." << endl;
						
			// Allocate memory in one big block for all square clusters
			int nNewClusters = (numRows-1)*(numCols-1);
			SquareCluster newCluster[nNewClusters];
			
			vector<int> sqAssignment; // temp to use below
			sqAssignment.push_back(-1); sqAssignment.push_back(-1); sqAssignment.push_back(-1); sqAssignment.push_back(-1);
			
			int index=0;
			for(int r=0; r < numRows-1; r++)
				for(int c=0; c < numCols-1; c++)
				{
					newCluster[index].i = sub2ind(numCols, numRows, c, r); // node which is in position (r,c)
					newCluster[index].j = sub2ind(numCols, numRows, c+1, r); // node which is in position (r,c+1)
					newCluster[index].k = sub2ind(numCols, numRows, c+1, r+1); // node which is in position (r+1,c+1)
					newCluster[index].l = sub2ind(numCols, numRows, c, r+1); // node which is in position (r+1,c)
					
					// Find the intersection sets for this square
					vector<int> ij_edge; ij_edge.push_back(newCluster[index].i); ij_edge.push_back(newCluster[index].j);
					newCluster[index].ij_intersect_loc = mplp.FindIntersectionSet(ij_edge);
					vector<int> kj_edge; kj_edge.push_back(newCluster[index].k); kj_edge.push_back(newCluster[index].j);
					newCluster[index].kj_intersect_loc = mplp.FindIntersectionSet(kj_edge);
					vector<int> kl_edge; kl_edge.push_back(newCluster[index].k); kl_edge.push_back(newCluster[index].l);
					newCluster[index].kl_intersect_loc = mplp.FindIntersectionSet(kl_edge);
					vector<int> il_edge; il_edge.push_back(newCluster[index].i); il_edge.push_back(newCluster[index].l);
					newCluster[index].il_intersect_loc = mplp.FindIntersectionSet(il_edge);
					
					vector<MulDimArr*> beliefs;
					vector<bool> b_transpose;
					beliefs.push_back(&mplp.m_sum_into_intersects[newCluster[index].ij_intersect_loc]);
					b_transpose.push_back(mplp.m_all_intersects[newCluster[index].ij_intersect_loc][0] != newCluster[index].i); // i first
					beliefs.push_back(&mplp.m_sum_into_intersects[newCluster[index].kj_intersect_loc]);
					b_transpose.push_back(mplp.m_all_intersects[newCluster[index].kj_intersect_loc][0] != newCluster[index].j); // then 'j'
					beliefs.push_back(&mplp.m_sum_into_intersects[newCluster[index].kl_intersect_loc]);
					b_transpose.push_back(mplp.m_all_intersects[newCluster[index].kl_intersect_loc][0] != newCluster[index].k); // then 'k'
					beliefs.push_back(&mplp.m_sum_into_intersects[newCluster[index].il_intersect_loc]);
					b_transpose.push_back(mplp.m_all_intersects[newCluster[index].il_intersect_loc][0] != newCluster[index].l); // then 'l'
					
					double bound_indep = maximizeIndependently(beliefs);
					
					// Before doing expensive joint maximization, see if we can quickly find optimal assignment
					sqAssignment[0] = mplp.m_decoded_res[newCluster[index].i]; sqAssignment[1] = mplp.m_decoded_res[newCluster[index].j];
					sqAssignment[2] = mplp.m_decoded_res[newCluster[index].k]; sqAssignment[3] = mplp.m_decoded_res[newCluster[index].l];
					double bound_quick = getValCycle(beliefs, b_transpose, sqAssignment);
					if(bound_indep == bound_quick)
					{
						newCluster[index].bound = 0;
					} else {
						// Do expensive joint maximization
						// Calculate bound for the square whose upper left corner is grid position (r,c)
						newCluster[index].bound = bound_indep - maximizeCycle(beliefs, b_transpose);
					}
					
					index++;
				}

			// TODO opt: have a class for a cluster, so we can have different types and sort by bound,
			//       choosing best one by bound.
			//       Make the sorting and adding independent of the type of graph...
			
			// Sort the clusters by the bound
			sort(newCluster, newCluster+nNewClusters);
			
			printf(" -- Considered %d clusters, smallest bound %g, largest bound %g\n", nNewClusters, newCluster[nNewClusters-nclus_to_add].bound, newCluster[nNewClusters-1].bound);
			
			// Add the top nclus_to_add clusters to the relaxation
			for(int clusterId = nNewClusters-1; clusterId >= 0 && nClustersAdded < nclus_to_add; clusterId--)
			{
				// TODO: check that these clusters and intersection sets haven't already been added

				// Triangulate along edge (i,k)

				// Need to add new intersection set
				vector<int> ik_edge; ik_edge.push_back(newCluster[clusterId].i); ik_edge.push_back(newCluster[clusterId].k);
				int ik_intersect_loc = mplp.AddIntersectionSet(ik_edge);

				// Now add clusters ijk and ilk
				vector<int> ijk_inds, ilk_inds;
				ijk_inds.push_back(newCluster[clusterId].i); ijk_inds.push_back(newCluster[clusterId].j); ijk_inds.push_back(newCluster[clusterId].k);
				ilk_inds.push_back(newCluster[clusterId].i); ilk_inds.push_back(newCluster[clusterId].l); ilk_inds.push_back(newCluster[clusterId].k);
				
				vector<int> ijk_intersect_inds, ilk_intersect_inds;
				ijk_intersect_inds.push_back(ik_intersect_loc);
				ijk_intersect_inds.push_back(newCluster[clusterId].ij_intersect_loc);
				ijk_intersect_inds.push_back(newCluster[clusterId].kj_intersect_loc);
				
				ilk_intersect_inds.push_back(ik_intersect_loc);
				ilk_intersect_inds.push_back(newCluster[clusterId].il_intersect_loc);
				ilk_intersect_inds.push_back(newCluster[clusterId].kl_intersect_loc);
				
				mplp.AddRegion(ijk_inds, ijk_intersect_inds);
				mplp.AddRegion(ilk_inds, ilk_intersect_inds);

				// TODO: log which clusters are chosen...
				
				nClustersAdded++;
			}
		}
		
		finish = clock();
		time = (double(finish)-double(start))/CLOCKS_PER_SEC;
		printf(" -- Added %d clusters to relaxation. Took %lg seconds\n", nClustersAdded, time);
		
		// TODO: log different algo_triplet iterations, as we did in Matlab. Note that obj_del is 1e10.
		
		// Run MPLP again
		printf("Running MPLP again for %d more iterations\n", niter_later);
		mplp.RunMPLP(niter_later,obj_del_thr,int_gap_thr);
	}
	
	// TODO opt: checkpointing... so we have this printed more often
	printf("Writing final results to output files...\n");
	mplp.Write("msgs.txt","res.txt","suminto.txt","objhist.txt","inthist.txt","timehist.txt");

	return 1;
}
