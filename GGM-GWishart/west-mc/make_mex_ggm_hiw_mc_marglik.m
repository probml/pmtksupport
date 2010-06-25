function make_mex_ggm_hiw_mc_marglik()

fprintf('make_mex_ggm_hiw_mc_marglik: mex file not found. Attempting to compile...')

if(strcmp(computer,'GLNX86') | strcmp(computer,'GLNXA64'))
  mex -g -L/usr/lib -largeArrayDims -lmwblas -lmwlapack  -lm mex_ggm_hiw_mc_marglik.cpp newgraph.cpp variousfunctions.cpp likecomb.cpp  ISampling.cpp random.cpp graphfns.cpp 
elseif(strcmp(computer,'PCWIN'))
  liblapack = fullfile(matlabroot,'extern', 'lib', 'win32', 'microsoft','libmwlapack.lib');
  libblas = fullfile(matlabroot,'extern', 'lib', 'win32', 'microsoft','libmwblas.lib');
  mex('-g','-largeArrayDims','mex_ggm_hiw_mc_marglik.cpp','newgraph.cpp','variousfunctions.cpp','likecomb.cpp','ISampling.cpp','random.cpp','graphfns.cpp',liblapack,libblas);
elseif(strcmp(computer,'PCWIN64'))
  liblapack = fullfile(matlabroot,'extern', 'lib', 'win64', 'microsoft','libmwlapack.lib');
  libblas = fullfile(matlabroot,'extern', 'lib', 'win64', 'microsoft','libmwblas.lib');
  mex('-g','-largeArrayDims','mex_ggm_hiw_mc_marglik.cpp','newgraph.cpp','variousfunctions.cpp','likecomb.cpp','ISampling.cpp','random.cpp','graphfns.cpp',liblapack,libblas);
else
  fprintf('  Warning: This code package has not been tested on your platform (%s).\n',computer);
  fprintf('  Default mex options and library paths may not be set correctly in make_mex_ggm_hiw_mc_marglik.\n');
  fprintf('  Please create a configuration for your platform in make_mex_ggm_hiw_mc_marglik.m.\n');
  edit('make_mex_ggm_hiw_mc_marglik.m');
end

