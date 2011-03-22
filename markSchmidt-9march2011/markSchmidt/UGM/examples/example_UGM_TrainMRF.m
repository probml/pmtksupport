%%
clear all
close all
load rain.mat
y = X+1; % Convert from {0,1} to {1,2} representation

% Plot what data looks like
figure(1);
imagesc(y(1:100,:))
title('Rain Data for first 100 months');

% Compute marginal of raining on any day
p_rain = sum(y(:)==2)/numel(y)

% Compute log-likelihood of full data set
negloglik_y = log(p_rain)*sum(y(:)==2) + log(1-p_rain)*sum(y(:)==1)

% Plot what independent samples would look like
figure(2);
imagesc(p_rain > rand(100,28));
title('Samples based on independent model');
fprintf('(paused)\n');
pause

%% Make edgeStruct
[nInstances,nNodes] = size(y);
nStates = max(y);
adj = zeros(nNodes);
for i = 1:nNodes-1
    adj(i,i+1) = 1;
end
adj = adj+adj';
edgeStruct = UGM_makeEdgeStruct(adj,nStates);
maxState = max(nStates);
nEdges = edgeStruct.nEdges;

%% Training (Ising = 1)

% Make nodeMap
nodeMap = zeros(nNodes,maxState);
nodeMap(:,1) = 1;

% Make edgeMap
edgeMap = zeros(maxState,maxState,nEdges);
edgeMap(1,1,:) = 2;
edgeMap(2,2,:) = 2;

% Initialize weights
nParams = max([nodeMap(:);edgeMap(:)]);
w = zeros(nParams,1);

% Example of making potentials
[nodePot,edgePot] = UGM_MRF_makePotentials(w,nodeMap,edgeMap,edgeStruct);

% Compute sufficient statistics
suffStat = UGM_MRF_computeSuffStat(y,nodeMap,edgeMap,edgeStruct);

% Evaluate NLL
nll = UGM_MRF_NLL(w,nInstances,suffStat,nodeMap,edgeMap,edgeStruct,@UGM_Infer_Chain)

% Optimize
w = minFunc(@UGM_MRF_NLL,w,[],nInstances,suffStat,nodeMap,edgeMap,edgeStruct,@UGM_Infer_Chain)

% Now make potentials
[nodePot,edgePot] = UGM_MRF_makePotentials(w,nodeMap,edgeMap,edgeStruct);
nodePot(1,:)
edgePot(:,:,1)
fprintf('(paused)\n');
pause

%% Training (Ising = 2)

edgeMap(2,2,:) = 3;

% Initialize weights
nParams = max([nodeMap(:);edgeMap(:)]);
w = zeros(nParams,1);

% Example of making potentials
[nodePot,edgePot] = UGM_MRF_makePotentials(w,nodeMap,edgeMap,edgeStruct);

% Compute sufficient statistics
suffStat = UGM_MRF_computeSuffStat(y,nodeMap,edgeMap,edgeStruct);

% Optimize
w = minFunc(@UGM_MRF_NLL,w,[],nInstances,suffStat,nodeMap,edgeMap,edgeStruct,@UGM_Infer_Chain)

% Now make potentials
[nodePot,edgePot] = UGM_MRF_makePotentials(w,nodeMap,edgeMap,edgeStruct);
nodePot(1,:)
edgePot(:,:,1)
fprintf('(paused)\n');
pause

%% Training (Ising = 0)

edgeMap(1,2,:) = 4;

% Initialize weights
nParams = max([nodeMap(:);edgeMap(:)]);
w = zeros(nParams,1);

% Example of making potentials
[nodePot,edgePot] = UGM_MRF_makePotentials(w,nodeMap,edgeMap,edgeStruct);

% Compute sufficient statistics
suffStat = UGM_MRF_computeSuffStat(y,nodeMap,edgeMap,edgeStruct);

% Optimize
w = minFunc(@UGM_MRF_NLL,w,[],nInstances,suffStat,nodeMap,edgeMap,edgeStruct,@UGM_Infer_Chain)

% Now make potentials
[nodePot,edgePot] = UGM_MRF_makePotentials(w,nodeMap,edgeMap,edgeStruct);
nodePot(1,:)
edgePot(:,:,1)
fprintf('(paused)\n');
pause

%% Do decoding/infence/sampling in learned model

decode = UGM_Decode_Chain(nodePot,edgePot,edgeStruct)

[nodeBel,edgeBel,logZ] = UGM_Infer_Chain(nodePot,edgePot,edgeStruct);
nodeBel

samples = UGM_Sample_Chain(nodePot,edgePot,edgeStruct);
figure(3);
imagesc(samples')
title('Samples from MRF model');
fprintf('(paused)\n');
pause

%% Do conditional decoding/inference/sampling in learned model

clamped = zeros(nNodes,1);
clamped(1:2) = 2;

condDecode = UGM_Decode_Conditional(nodePot,edgePot,edgeStruct,clamped,@UGM_Decode_Chain)
condNodeBel = UGM_Infer_Conditional(nodePot,edgePot,edgeStruct,clamped,@UGM_Infer_Chain)
condSamples = UGM_Sample_Conditional(nodePot,edgePot,edgeStruct,clamped,@UGM_Sample_Chain);

figure(4);
imagesc(condSamples')
title('Conditional samples from MRF model');
fprintf('(paused)\n');
pause