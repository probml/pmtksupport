/*
 *  mplp_main.cpp
 *  mplp
 *
 *  Created by Amir Globerson on 6/10/08.
 *  Copyright 2008 MIT. All rights reserved.
 *
 */
#include <iostream>

using namespace std;
#include <vector>
#include "mplp_alg.h"

int main( int argc, char *argv[] )
{
	int nIter = 10;
	double obj_del_thr;
	double int_gap_thr;
	
	if (argc<4)
	{
		printf("Syntax: cmplp niter obj_del_thr int_gap_thr\n");
		return 0;
	}
	
	sscanf(argv[1],"%d",&nIter);
	sscanf(argv[2],"%lg",&obj_del_thr);
	sscanf(argv[3],"%lg",&int_gap_thr);	
	
	printf("nIters=%d\nobj_del_thr=%lg\nint_gap_thr=%lg\n",nIter,obj_del_thr,int_gap_thr);

	MPLPAlg mplp("regions.txt","intersects.txt","region_intersects.txt","lambdas.txt","var_sizes.txt");
	mplp.Read("msgs.txt");
		
	mplp.RunMPLP(nIter,obj_del_thr,int_gap_thr);
	
	// We always write the output
	mplp.Write("msgs.txt","res.txt","suminto.txt","objhist.txt","inthist.txt","timehist.txt");
	return 1;
}
