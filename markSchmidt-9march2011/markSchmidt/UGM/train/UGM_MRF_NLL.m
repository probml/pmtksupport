function [NLL,g] = UGM_MRF_NLL(w,nInstances,suffStat,nodeMap,edgeMap,edgeStruct,inferFunc,varargin)

[nNodes,maxState] = size(nodeMap);
nEdges = edgeStruct.nEdges;
edgeEnds = edgeStruct.edgeEnds;
nStates = edgeStruct.nStates;

% Make potentials
[nodePot,edgePot] = UGM_MRF_makePotentials(w,nodeMap,edgeMap,edgeStruct);

% Compute marginals and logZ
[nodeBel,edgeBel,logZ] = inferFunc(nodePot,edgePot,edgeStruct,varargin{:});

% Compute NLL
NLL = -w'*suffStat + nInstances*logZ;

if nargout > 1
    g = -suffStat;
    for n = 1:nNodes
        for s = 1:nStates(n)
            if nodeMap(n,s) > 0
                g(nodeMap(n,s)) = g(nodeMap(n,s)) + nInstances*nodeBel(n,s);
            end
        end
    end
    for e = 1:nEdges
        n1 = edgeEnds(e,1);
        n2 = edgeEnds(e,2);
        for s1 = 1:nStates(n1)
            for s2 = 1:nStates(n2)
                if edgeMap(s1,s2,e) > 0
                    g(edgeMap(s1,s2,e)) = g(edgeMap(s1,s2,e)) + nInstances*edgeBel(s1,s2,e);
                end
            end
        end
    end
end