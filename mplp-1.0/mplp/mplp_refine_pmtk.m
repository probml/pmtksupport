% Modifications to mplp_refine for the pmtk package

% Function: mplp_refine
% Implements a MAP-LP approximation algorithm as described in
% Sontag et al, UAI 08.
% 
% The input defines a graphical model over n variables.
%
% Input:
%        adj    - An adjacency matrix for a graph, where adj(i,j)=1 if
%                 there is an edge between i and j, and zero otherwise.
%
%        lambda - The interaction potentials for the graphical
%                 model. lambda{xi,xj} is matrix m where m(k,l) is
%                 theta_{ij}(xi=k,xj=l) (see notations in the paper)
%
%        local  - The local potentials for the graphical
%                 model. local{xi} is a vector v where v(k) is theta_i(xi) 
%
%            [[ Taken together these define a function f(x) =
%               \sum_{ij:adj(i,j)=1}\theta_{ij}(xi,xj)+\sum_i\theta_i(xi)
%             ]]
%       
%        params - A structure with the following fields (default
%                 values in brackets):
%
%                 file_type (0): Determines the regions that the
%                          algorithm tries to add to the approximation.
%                          If file_type=0 it adds triplets regions on
%                          triangles that already exist in the graph.
%                          If file_type=1 it assumes the graph is a grid,
%                          and adds square regions (see more below).
%                   
%                 niter (1000): The number of iterations to run
%                               MPLP on the model with only edge regions 
%                 niter_later (20): The number of iterations to run
%                               MPLP after each batch of added regions  
%                 nclust_to_add (5): The number of regions to add
%                                    in each batch 
%                 obj_del_thr (1e-4): The algorithm will stop if
%                      the objective changes by less than obj_del_thr in
%                      two consecutive iterations      
%                 int_gap_thr (1e-5): The algorithm will stop if
%                      the gap between dual and primal objective is less
%                      than int_gap_thr. In this case it will return the
%                      exact MAP.
%                 b_do_max (1): If b_do_max=1 the objective f(x) is
%                      maximized, if b_do_max=0 it is minimized
%
%                 The following two parameters are only used for
%                 grid graphs. NOTE: for grid graphs, the nodes
%                 must be input in row-major order. (This is so that
%                 we know where to look for the squares, when tightening
%                 the relaxation.)
%
%                 num_rows (sqrt(#nodes)): For grid graphs, the number
%                                          of rows in the graph.
%                 num_cols (sqrt(#nodes)): For grid graphs, the number
%                                          of rows in the graph.
%
% Output:
%       assign - The optimal assignment found by the algorithm
%       assign_val - The model value for assign (i.e., f(assign)) 
%       dual_obj_hist - Values of the dual objective for each iteration
%       int_val_hist  - Values of the model "probability" for the
%                       integral assignments decoded at each iteration.  
%       NOTE: dual_obj_hist is always greater than int_val_hist. If
%       these values coincide, the exact MAP was found (this is for
%       maximization problems, for minimization it is the opposite).

function [assign,dual_obj_hist,int_val_hist] = mplp_refine(adj,lambda,local,params)
% Set default values for params
if ~exist('params','var')
  params = [];
end
% Will set the values for fields not already in params
params = set_default_params(params,{'file_type','niter', ...
                    'niter_later','nclust_to_add','obj_del_thr', ...
                    'int_gap_thr','b_do_max'},[0 1000 20 5 1e-4 1e-5 1],length(local));


%params

% Negate field for minimization
if ~params.b_do_max
   [lambda,local] = negate_field(adj,lambda,local);
end

[row_inds,col_inds] = find(triu(adj));
numPairs = length(row_inds);
numNodes = size(adj,1);

regions = cell(1,numPairs+numNodes);
regions(1:numPairs) = mat2cell([row_inds(:), col_inds(:)],ones(1,numPairs),2);
regions(numPairs+1:numPairs+numNodes) = mat2cell([1:numNodes]', ones(1,numNodes), 1);

% We need the edge intersection sets to calculate the bound criterion
intersects = regions;

my_region_lambda = {};

for i=1:numPairs
  tmp = lambda{row_inds(i),col_inds(i)}';
  my_region_lambda{i} = tmp(:)';
end
for i=1:numNodes
  my_region_lambda{end+1} = local{i}';
end

% We no longer need lambda, so set to 0.
lambda = 0;

      
% Initialize MPLP with these regions and intersections
gmplp_state = gmplp_init(regions,intersects,my_region_lambda,local,-1,0,lambda);

time_init = toc;


% we save the temporary data to the directory
% where the executable is so it can find it
folder = fullfile(pmtk3Root(), 'data')
%[folder, fname, extension] = fileparts(which('mplp_refine_pmtk')); %#ok
fprintf('Saving gmplp_state to %s\n', folder);


save_gmplp_state_for_c(gmplp_state, folder);

fprintf('Done Saving\n');

if params.file_type == 0
  
  % Tighten with triplets (e.g. for protein files)
  
  % Call c code
  str = sprintf('algo_triplet %d %d %d %g %g 0',params.niter, ...
    params.niter_later,params.nclust_to_add, ...
    params.obj_del_thr,params.int_gap_thr);
  
elseif params.file_type == 1
  
  % Tighten with squares (grids)
  
  % Call c code
  str= sprintf('algo_triplet %d %d %d %g %g 1 %d %d',params.niter, ...
    params.niter_later,params.nclust_to_add, ...
    params.obj_del_thr,params.int_gap_thr, ...
    params.num_rows, params.num_cols);
end


if 0
  eval(sprintf('!./%s', str))
else
  % ensure that any temporary files are written 
  % to the location of the executable so it can find them
  curFolder = pwd;
  cd(folder)
  system(str)
  cd(curFolder)
end


load(fullfile(folder, 'inthist.txt'))
load(fullfile(folder, 'objhist.txt'))
load(fullfile(folder, 'res.txt'))

assign = res+1;
dual_obj_hist = objhist;
int_val_hist = inthist;

if ~params.b_do_max
  dual_obj_hist = -dual_obj_hist;
  int_val_hist = -int_val_hist;
end

assign_val = int_val_hist(end);

end
