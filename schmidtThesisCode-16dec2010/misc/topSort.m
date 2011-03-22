function [Q] = topSort(adj)

nNodes = length(adj);
Q = zeros(nNodes,1);
endQ = 1;
for i = 1:nNodes
    if all(adj(:,i)==0)
        Q(endQ) = i;
        endQ = endQ + 1;
    end
end
nEdges = sum(adj);
for i = 1:nNodes
    for j = 1:nNodes
        if adj(Q(i),j) == 1
            adj(Q(i),j) = 0;
            nEdges(j) = nEdges(j)-1;
            if nEdges(j) == 0
                Q(endQ) = j;
                endQ = endQ+1;
            end
        end
    end
end

if 1 % Check correctness
   for i = 1:nNodes
       if any(adj(Q(1:i-1),i)==1)
           fprintf('Topsort did not work\n');
           pause
       end
   end
   if length(unique(Q)) ~= nNodes
       fprintf('Topsort did not work\n');
   end
end