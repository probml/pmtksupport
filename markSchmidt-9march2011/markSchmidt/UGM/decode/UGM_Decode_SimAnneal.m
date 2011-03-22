function  [y] = UGM_Infer_Gibbs(nodePot, edgePot, edgeStruct,y)
% INPUT
% nodePot(node,class)
% edgePot(class,class,edge) where e is referenced by V,E (must be the same
% between feature engine and inference engine)
%
% OUTPUT
% nodeLabel(node)

[nNodes,maxStates] = size(nodePot);
nEdges = size(edgePot,3);
edgeEnds = edgeStruct.edgeEnds;
V = edgeStruct.V;
E = edgeStruct.E;
nStates = edgeStruct.nStates;
maxIter = edgeStruct.maxIter;

if nargin < 4
    % Initialize
    [junk y] = max(nodePot,[],2);
end

% Compute rate of temperature decrease
maxTemp = 10;
minTemp = 1e-2;
Trate = exp((log(minTemp)-log(maxTemp))/maxIter);

T = maxTemp;
for i = 1:maxIter
    for n = 1:nNodes
        y_old = y(n);

        % Generate Proposal (uniform over other possible states)
        s = ceil(rand*(nStates(n)-1));
        if s >= y(n)
            y(n) = s+1;
        else
            y(n) = s;
        end

        % Compute Node Potential
        pot = log(nodePot(n,[y_old y(n)]));

        % Find Neighbors
        edges = E(V(n):V(n+1)-1);

        % Multiply Edge Potentials
        for e = edges(:)'
            n1 = edgeEnds(e,1);
            n2 = edgeEnds(e,2);

            if n == edgeEnds(e,1)
                ep = edgePot([y_old y(n)],y(n2),e)';
            else
                ep = edgePot(y(n1),[y_old y(n)],e);
            end
            pot = pot + log(ep);
        end
        
        % Compute Metropolis ratio
        metRatio = exp((pot(2)-pot(1))/T);

        if rand > metRatio
            y(n) = y_old; % Reject move
        end
    end
    % Update Temperature
    T = T*Trate;
end