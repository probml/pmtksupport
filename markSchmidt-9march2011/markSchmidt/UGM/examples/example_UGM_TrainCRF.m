%%
clear all
close all
load rain.mat

% Make rain labels y, and binary month features X
y = X+1;
[nInstances,nNodes] = size(y);

%% Make edgeStruct
nStates = max(y);
adj = zeros(nNodes);
for i = 1:nNodes-1
    adj(i,i+1) = 1;
end
adj = adj+adj';
edgeStruct = UGM_makeEdgeStruct(adj,nStates);
nEdges = edgeStruct.nEdges;
maxState = max(nStates);

%% Temp
edgeStruct.nStates = int32(edgeStruct.nStates);
edgeStruct.edgeEnds = int32(edgeStruct.edgeEnds);
y = int32(y);

%% Training (no features)

% Make simple bias features
Xnode = ones(nInstances,1,nNodes);
Xedge = ones(nInstances,1,nEdges);

% Make nodeMap
nodeMap = zeros(nNodes,maxState,'int32');
nodeMap(:,1) = 1;

edgeMap = zeros(maxState,maxState,nEdges,'int32');
edgeMap(1,1,:) = 2;
edgeMap(2,1,:) = 3;
edgeMap(1,2,:) = 4;

% Initialize weights
nParams = max([nodeMap(:);edgeMap(:)]);
w = zeros(nParams,1);

% Example of making potentials for the first training example
instance = 1;
[nodePot,edgePot] = UGM_CRF_makePotentials(w,Xnode,Xedge,nodeMap,edgeMap,edgeStruct,instance);

% Optimize
w = minFunc(@UGM_CRF_NLL,randn(size(w)),[],Xnode,Xedge,y,nodeMap,edgeMap,edgeStruct,@UGM_Infer_Chain)
fprintf('(paused)\n');
pause

%% Training (with node features, but no edge features)

% Make simple bias features
nFeatures = 12;
Xnode = zeros(nInstances,nFeatures,nNodes);
for m = 1:nFeatures
    Xnode(months==m,m,:) = 1;
end
Xnode = [ones(nInstances,1,nNodes) Xnode];
nNodeFeatures = size(Xnode,2);

% Make nodeMap
nodeMap = zeros(nNodes,maxState,nNodeFeatures,'int32');
for f = 1:nNodeFeatures
    nodeMap(:,1,f) = f;
end

% Make edgeMap
edgeMap = zeros(maxState,maxState,nEdges,'int32');
edgeMap(1,1,:) = nFeatures+1;
edgeMap(2,1,:) = nFeatures+2;
edgeMap(1,2,:) = nFeatures+3;

% Initialize weights
nParams = max([nodeMap(:);edgeMap(:)]);
w = zeros(nParams,1);

% Optimize
w = minFunc(@UGM_CRF_NLL,w,[],Xnode,Xedge,y,nodeMap,edgeMap,edgeStruct,@UGM_Infer_Chain)
fprintf('(paused)\n');
pause

%% Training (with edge features)

% Make edge features
sharedFeatures = 1:13;
Xedge = UGM_makeEdgeFeatures(Xnode,edgeStruct.edgeEnds,sharedFeatures);
nEdgeFeatures = size(Xedge,2);

% Make nodeMap
nodeMap = zeros(nNodes,maxState,nNodeFeatures,'int32');
for f = 1:nNodeFeatures
    nodeMap(:,1,f) = f;
end

% Make edgeMap
for edgeFeat = 1:nEdgeFeatures
    for s1 = 1:2
        for s2 = 1:2
            f = f+1;
            edgeMap(s1,s2,:,edgeFeat) = f;
        end
    end
end

% Initialize weights
nParams = max([nodeMap(:);edgeMap(:)]);
w = zeros(nParams,1);

% Optimize
UGM_CRF_NLL(w,Xnode,Xedge,y,nodeMap,edgeMap,edgeStruct,@UGM_Infer_Chain);
w = minFunc(@UGM_CRF_NLL,w,[],Xnode,Xedge,y,nodeMap,edgeMap,edgeStruct,@UGM_Infer_Chain)
fprintf('(paused)\n');
pause

%% Do decoding/infence/sampling in learned model (given features)

% We will look at a case in December
i = 11;
[nodePot,edgePot] = UGM_CRF_makePotentials(w,Xnode,Xedge,nodeMap,edgeMap,edgeStruct,i);

decode = UGM_Decode_Chain(nodePot,edgePot,edgeStruct)

[nodeBel,edgeBel,logZ] = UGM_Infer_Chain(nodePot,edgePot,edgeStruct);
nodeBel

samples = UGM_Sample_Chain(nodePot,edgePot,edgeStruct);
figure(1);
imagesc(samples')
title('Samples from CRF model (for December)');
fprintf('(paused)\n');
pause

%% Do conditional decoding/inference/sampling in learned model (given features)

clamped = zeros(nNodes,1);
clamped(1:2) = 2;

condDecode = UGM_Decode_Conditional(nodePot,edgePot,edgeStruct,clamped,@UGM_Decode_Chain)
condNodeBel = UGM_Infer_Conditional(nodePot,edgePot,edgeStruct,clamped,@UGM_Infer_Chain)
condSamples = UGM_Sample_Conditional(nodePot,edgePot,edgeStruct,clamped,@UGM_Sample_Chain);

figure(2);
imagesc(condSamples')
title('Conditional samples from CRF model (for December)');
fprintf('(paused)\n');
pause

%% Now see what samples in July look like

XtestNode = [1 0 0 0 0 0 0 1 0 0 0 0 0]; % Turn on bias and indicator variable for July
XtestNode = repmat(XtestNode,[1 1 nNodes]);
XtestEdge = UGM_makeEdgeFeatures(XtestNode,edgeStruct.edgeEnds,sharedFeatures);

[nodePot,edgePot] = UGM_CRF_makePotentials(w,XtestNode,XtestEdge,nodeMap,edgeMap,edgeStruct);

samples = UGM_Sample_Chain(nodePot,edgePot,edgeStruct);
figure(3);
imagesc(samples')
title('Samples from CRF model (for July)');
fprintf('(paused)\n');
pause

%% Training with L2-regularization

% Set up regularization parameters
lambda = 10*ones(size(w));
lambda(1) = 0; % Don't penalize node bias variable
lambda(14:17) = 0; % Don't penalize edge bias variable
regFunObj = @(w)penalizedL2(w,@UGM_CRF_NLL,lambda,Xnode,Xedge,y,nodeMap,edgeMap,edgeStruct,@UGM_Infer_Chain);

% Optimize
w = zeros(nParams,1);
w = minFunc(regFunObj,w);
UGM_CRF_NLL(w,Xnode,Xedge,y,nodeMap,edgeMap,edgeStruct,@UGM_Infer_Chain);