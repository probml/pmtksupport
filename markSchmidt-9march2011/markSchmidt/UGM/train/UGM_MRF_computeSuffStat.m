function [suffStat] = UGM_MRF_computeSuffStat(Y,nodeMap,edgeMap,edgeStruct)

[nNodes,maxState] = size(nodeMap);
nEdges = edgeStruct.nEdges;
edgeEnds = edgeStruct.edgeEnds;
nStates = edgeStruct.nStates;
nParams = max([nodeMap(:);edgeMap(:)]);

nInstances = size(Y,1);
suffStat = zeros(nParams,1);
for i = 1:nInstances
   y = Y(i,:);
   for n = 1:nNodes
      if nodeMap(n,y(n)) > 0
         suffStat(nodeMap(n,y(n))) = suffStat(nodeMap(n,y(n))) + 1;
      end
   end
   for e = 1:nEdges
      n1 = edgeEnds(e,1);
      n2 = edgeEnds(e,2);
      if edgeMap(y(n1),y(n2),e) > 0
         suffStat(edgeMap(y(n1),y(n2),e)) = suffStat(edgeMap(y(n1),y(n2),e)) + 1;
      end
   end
end