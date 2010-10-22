niter = 1000;       % Number of MPLP iterations initially.
niter_later = 100;  % Number of iters after adding each square.
nclust_to_add = 5;  % Add this many squares before re-solving.
obj_del_thr = 2e-4; % Stopping criterion.
int_gap_thr = 2e-4;

% Choose some parameters for generating the random grid
seed = 0;
coupling = 2;
field = 1;
nPerSide = 10;
nVals = 3;

% Fix random seed for reproducibility
rand('twister', seed)

% Generate the random grid (the test case)
% NOTE: this is a bit inefficient, often slower than finding MAP :)
[G,local,lambda,exp_local,exp_lambda,xs,s_pair,Fields] = tst_potts('grid',nVals,nPerSide,3,coupling,field);
adj = G.Edges;

numNodes=size(adj,1);

% Call MPLP

params = struct;
params.file_type = 1; % for grid
params.niter = niter;
params.niter_later = niter_later;
params.nclust_to_add = nclust_to_add;
params.obj_del_thr = obj_del_thr;
params.int_gap_thr = int_gap_thr;

% These are set by default to these values:
%params.num_cols = sqrt(numNodes);
%params.num_rows = sqrt(numNodes);

[assign,dual_obj_hist,int_val_hist] = mplp_refine(adj,lambda,local,params);




% Show the MAP assignment
assign

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
