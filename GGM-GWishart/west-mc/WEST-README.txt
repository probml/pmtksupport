This directory contains a program that uses a metropolis hastings
algorithm to sample general graphs according to their (approximate) posterior
probability. This file describes the compilation of the program, the 
input parameters to be set by the
user in an input file, the output format, and other parameters
advanced users can modify by changing and recompiling the source code.

COMPILATION: 
Set the appropriate library paths at line 7 of the Makefile.
Type "make metropolis".

INPUT:
The following parameters are specified by the user and read in via the 
standard input (eg ./metropolis < infile ).  See the "infile" given
for an example.
1) Data File (character string, no whitespace, not to exceed 512 characters)
2) Graph File (character string, no whitespace, not to exceed 512 characters)
3) number of variables (columns in Data File)
4) number of observations (rows in Data File)
5) tau
6) delta
7) number of metropolis Hastings iterations
8) random number seed

5) and 6) specify the parameters of the local 
inverse-Wishart prior over parameters used to obtain the graph
marginal likelihoods (prior described in Roverato 2002).  
The inverse-Wishart has degrees of freedom delta 
(should be >2 to ensure non-degeneracy) and matrix parameter tau*I. 

The data in the data file should be centered before input. Because the Wishart
parameter is of the form tau*I, the data should also be standardized
or believed to be on a common scale.  An appropriate tau must be
selected; if tau is too large it will dwarf the effect of the data. 
The marginal prior mode for the variance of each variable  is tau/(delta+1);
we use this quantity to set an appropriate value for tau. For example, 
if the data has been standardized so all the variances are 1.0, tau 
would be set to delta+1.

The Graph File specifies a starting graph. Examples are in "empty.dat"
and "nonempty.dat". The Graph file starts with the number of lines in the file.
Each subsequent line specifies and edge or an isolated node. An edge
is preceeded by a 2 and then lists the edege endpoints.  An isolated
node is preceeded by a 1. 

OUTPUT
For each metropolis Hastings iteration, the estimated posterior
probability is printed, and the current graph.  The graph is specified
as a sequence of 0's and 1's. These are the strict upper triangle (no
diagonal) of the incidence matrix, listed row by row. 

OTHER PARAMETERS
A sparsity encouraging prior over graphs is used; the prior probability
of a graph is a function of the number of edges in the graph and a
penalty parameter beta:  beta^(#edges)(1-beta)^(total possible edges -
#edges).  We have found fixing beta at 2/(Number of Nodes -1) to be
useful; it is fixed at this value within the program.  (Specified at
line 38 of metropolis.cpp)

If the graph is decomposable, the posterior is exact.  If the graph is
not decomposable, the marginal likelihood is estimated via the
algorithm of Atay-Kayis and Massam 2003.  We used (component size)^3/2 as
the number of iterations for estimating the prior normalising constant
and (component size)^3 *3/2 for the posterior normalising constant,
but with a minimum of  1000 iterations. This is specified at lines
139-147 of likecomb.cpp.  The minimum number of iterations should not
be changed; the first 1000 iterations, in addition to being used in
the normalizing constant estimation,  are used to guage the level of
the normalizing constant and prevent numerical problems with  the
exponentiation at line 246 of ISampling.cpp.
 
