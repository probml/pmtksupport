3/25/2009

This code is written by Amir Globerson and David Sontag, and
implements the algorithms described in the following two papers:

  Fixing max-product: Convergent message passing algorithms for MAP LP-relaxations
  Amir Globerson, Tommi Jaakkola
  Advances in Neural Information Processing Systems (NIPS) 21. Vancouver, Canada. 2007. 

  Tightening LP Relaxations for MAP using Message Passing
  David Sontag, Talya Meltzer, Amir Globerson, Tommi Jaakkola and Yair Weiss
  Uncertainty in Artificial Intelligence (UAI). Helsinki, Finland. 2008. 

MPLP is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by
the Free Software Foundation


MPLP is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
General Public License for more details.

You should have received a copy of the GNU General Public License
along with MPLP; see the file gpl.txt.  If not, write to the Free
Software Foundation, 59 Temple Place - Suite 330, Boston, MA
02111-1307, USA.

% --------------------------------------------------------------------
INSTALLATION

After you unpack mplp_ver1.tgz, you will have 3 directories:

   cmplp: Contains the C++ code for our algorithm.
   mplp:  Matlab wrapper, and also two test cases (see below).
   packages/@sparse_cell: Chen Yanover's sparse cell Matlab code,
                          used only by the protein test cases.

To begin, do the following:

1) Call "install_mplp.sh". This will compile the C++ code and
   put a binary called "algo_triplet" in the "mplp" directory.

2) Open Matlab and go into the "mplp" directory. Run "test_grid_c".
   This test case will create a 10x10 Ising grid and will find its
   MAP assignment.

3) Look at "mplp/mplp_refine.m", which is the primary interface for
   calling our algorithm. Look at "test_grid_c.m" and "test_protein_c.m"
   for two examples of how to call it.

It is also possible to run the C++ code without using the Matlab
wrapper, however we recommend you look at the Matlab files to see
what the required formats are for inputting the problem to our
algorithm.

% --------------------------------------------------------------------
SEARCHING FOR CLUSTERS

Our current implementation can tighten the relaxation by adding either
clusters on triplets of 3 variables, for general graphs, or by adding
clusters on the squares of 4 variables in grid graphs. For the former,
we only add clusters for 3 variables when all 3 of their edges already
exist in the MRF.

% --------------------------------------------------------------------
PROTEIN DESIGN

If you would like to try to run experiments on the protein side-chain
or protein design problems as reported in our UAI '08 paper, you can
download the following datasets:

   Linear Programming Relaxations and Belief Propagation -- An Empirical Study
   Chen Yanover, Talya Meltzer, Yair Weiss
   Journal of Machine Learning Research; 7(Sep):1887--1907, 2006.
   http://jmlr.csail.mit.edu/papers/volume7/yanover06a/data.html

  (301 MB) http://jmlr.csail.mit.edu/papers/volume7/yanover06a/Rosetta_SCP_Dataset.tgz
  (2.5 GB) http://jmlr.csail.mit.edu/papers/volume7/yanover06a/Rosetta_Design_Dataset.tgz

Then, edit "test_protein_c.m" to point to the relevant directories, and
call it!