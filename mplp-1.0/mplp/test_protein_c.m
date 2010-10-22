niter = 1000;       % Number of MPLP iterations initially.
niter_later = 20;   % Number of iters after adding each triplet cluster.
nclust_to_add = 5;  % Add this many triplet clusters before re-solving.
obj_del_thr = 2e-4; % Stopping criterion.
int_gap_thr = 2e-4;

% TODO TODO: Edit below 'input_dir' variables to point to right location.

% It's essential that we have sparse_cell before loading protein files
SPARSE_CELL_DIR = '../packages/'; % Chen Yanover's @sparse_cell class
addpath(SPARSE_CELL_DIR);



% Choose the input data problem
% -------------------------------------------------------------------
if 1 % Test SIDE-CHAIN
    
  input_dir = 'Rosetta/'; % Location of sidechain placement MRFs
    
  % Side-chain problems have '.mat' suffix
  file_suffix = '.mat';

  % Easier side-chain (small):
  %  prot_name = '1neu';
  % prot_name = '1qb7';

  % Hardest side-chain (big):
  prot_name = '1a8i';

else % Test DESIGN

  input_dir = 'rosetta_design/'; % Location of protein design MRFs

  % Design problems have '.dee.mat' suffixes
  file_suffix = '.dee.mat';

  % Smallest design:
  %prot_name = '1bx7';

  % Second smallest design (TRW can't solve):
  prot_name = '1en2';

  % ~1 hour test case (TRW can't solve):
  %prot_name = '1nkd';
  %prot_name='1igd';
%  prot_name = '1vfy';
  
  % Mid-sized design (TRW can't solve):
  % prot_name = '1psr';
  % prot_name = '1a8o';

  % Largest design (TRW can't solve):
  %prot_name = '1fpo';

end
% -------------------------------------------------------------------

input_file = sprintf('%s%s%s', input_dir, prot_name, file_suffix)
load(input_file);

adj = adjMatrix;
lambda = Eij;
local = Ei;

params = struct;
params.b_do_max = 0; % Protein file, so _minimize_ energy.
params.file_type = 0; % for default (tighten with triplets)
params.niter = niter;
params.niter_later = niter_later;
params.nclust_to_add = nclust_to_add;
params.obj_del_thr = obj_del_thr;
params.int_gap_thr = int_gap_thr;

[assign,dual_obj_hist,int_val_hist] = mplp_refine(adj,lambda,local,params);


% Plot the dual and primal objective values
leg = {};
plot(dual_obj_hist, 'Color',[0 0 1]);
leg{1} = 'Dual';
hold on
plot(int_val_hist, 'Color',[1 0 0]);
leg{2} = 'Primal (value of assignment)';
legend(leg, 'Location', 'NorthEast');
xlabel('MPLP iterations');
ylabel('Objective');
hold off
