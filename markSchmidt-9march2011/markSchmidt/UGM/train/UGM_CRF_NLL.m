function [NLL,g] = UGM_CRF_NLL(w,Xnode,Xedge,Y,nodeMap,edgeMap,edgeStruct,inferFunc,varargin)

[nNodes,maxState] = size(nodeMap);
nNodeFeatures = size(Xnode,2);
nEdgeFeatures = size(Xedge,2);
nEdges = edgeStruct.nEdges;
edgeEnds = edgeStruct.edgeEnds;
nStates = edgeStruct.nStates;

nInstances = size(Y,1);
NLL = 0;
g = zeros(size(w));

for i = 1:nInstances
    
    % Make potentials
    if edgeStruct.useMex
        [nodePot,edgePot] = UGM_CRF_makePotentialsC(w,Xnode,Xedge,nodeMap,edgeMap,nStates,edgeEnds,i);
    else
        [nodePot,edgePot] = UGM_CRF_makePotentials(w,Xnode,Xedge,nodeMap,edgeMap,edgeStruct,i);
    end
    
    % Compute marginals and logZ
    tmp = edgeStruct.useMex;
    edgeStruct.useMex = 0;
    [nodeBel,edgeBel,logZ] = inferFunc(nodePot,edgePot,edgeStruct,varargin{:});
    edgeStruct.useMex = tmp;
    
    % Update NLL
    if edgeStruct.useMex
        NLL = NLL - UGM_LogConfigurationPotentialC(Y(i,:),nodePot,edgePot,edgeEnds) + logZ;
        
        % Updates in-place
        UGM_CRF_NLLC(g,int32(i),nodeBel,edgeBel,edgeEnds,nStates,nodeMap,edgeMap,Xnode,Xedge,Y);
    else
        NLL = NLL - UGM_LogConfigurationPotential(Y(i,:),nodePot,edgePot,edgeEnds) + logZ;
        
        if nargout > 1
            for n = 1:nNodes
                for s = 1:nStates(n)
                    for f = 1:nNodeFeatures
                        if nodeMap(n,s,f) > 0
                            if s == Y(i,n)
                                obs = 1;
                            else
                                obs = 0;
                            end
                            g(nodeMap(n,s,f)) = g(nodeMap(n,s,f)) + Xnode(i,f,n)*(nodeBel(n,s) - obs);
                        end
                    end
                end
            end
            for e = 1:nEdges
                n1 = edgeEnds(e,1);
                n2 = edgeEnds(e,2);
                for s1 = 1:nStates(n1)
                    for s2 = 1:nStates(n2)
                        for f = 1:nEdgeFeatures
                            if edgeMap(s1,s2,e,f) > 0
                                if s1 == Y(i,n1) && s2 == Y(i,n2)
                                    obs = 1;
                                else
                                    obs = 0;
                                end
                                g(edgeMap(s1,s2,e,f)) = g(edgeMap(s1,s2,e,f)) + Xedge(i,f,e)*(edgeBel(s1,s2,e) - obs);
                            end
                        end
                    end
                end
            end
        end
    end
end