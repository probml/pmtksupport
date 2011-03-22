function [nodePot,edgePot] = UGM_CRF_makePotentials(w,Xnode,Xedge,nodeMap,edgeMap,edgeStruct,i)

if nargin < 7
    i = 1;
end

if edgeStruct.useMex
    [nodePot,edgePot] = UGM_CRF_makePotentialsC(w,Xnode,Xedge,nodeMap,edgeMap,edgeStruct.nStates,edgeStruct.edgeEnds,i);
else
    nNodes = size(nodeMap,1);
    maxState = size(nodeMap,2);
    nNodeFeatures = size(Xnode,2);
    nEdgeFeatures = size(Xedge,2);
    nEdges = edgeStruct.nEdges;
    edgeEnds = edgeStruct.edgeEnds;
    nStates = edgeStruct.nStates;
    
    nodePot = zeros(nNodes,maxState);
    for n = 1:nNodes
        for s = 1:nStates(n)
            for f = 1:nNodeFeatures
                if nodeMap(n,s,f) > 0
                    nodePot(n,s) = nodePot(n,s) + w(nodeMap(n,s,f))*Xnode(i,f,n);
                end
            end
            nodePot(n,s) = exp(nodePot(n,s));
        end
    end
    
    if nargout > 1
        edgePot = zeros(maxState,maxState,nEdges);
        for e = 1:nEdges
            n1 = edgeEnds(e,1);
            n2 = edgeEnds(e,2);
            for s1 = 1:nStates(n1)
                for s2 = 1:nStates(n2)
                    for f = 1:nEdgeFeatures
                        if edgeMap(s1,s2,e,f) > 0
                            edgePot(s1,s2,e) = edgePot(s1,s2,e) + w(edgeMap(s1,s2,e,f))*Xedge(i,f,e);
                        end
                    end
                    edgePot(s1,s2,e) = exp(edgePot(s1,s2,e));
                end
            end
        end
    end
end