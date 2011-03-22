function [nodePot,edgePot] = UGM_MRF_makePotentials(w,nodeMap,edgeMap,edgeStruct)

[nNodes,maxState] = size(nodeMap);
nEdges = edgeStruct.nEdges;
edgeEnds = edgeStruct.edgeEnds;
nStates = edgeStruct.nStates;

nodePot = zeros(nNodes,maxState);
for n = 1:nNodes
    for s = 1:nStates(n)
        if nodeMap(n,s) == 0
            nodePot(n,s) = 1;
        else
            nodePot(n,s) = exp(w(nodeMap(n,s)));
        end
    end
end

if nargout > 1
    edgePot = zeros(maxState,maxState,nEdges);
    for e = 1:nEdges
        n1 = edgeEnds(e,1);
        n2 = edgeEnds(e,2);
        for s1 = 1:nStates(n1)
            for s2 = 1:nStates(n2)
                if edgeMap(s1,s2,e) == 0
                    edgePot(s1,s2,e) = 1;
                else
                    edgePot(s1,s2,e) = exp(w(edgeMap(s1,s2,e)));
                end
            end
        end
    end
end