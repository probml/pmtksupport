function [G,x,s_pair] = make_grid_graph(nPerSideX,nPerSideY,interact_sign,interact_strength,field_strength)


%nPerSide = 2;%
%edges = rand(nVertices,nVertices)>0.5;
x=[];
%rand('seed',0);
for i=1:nPerSideX
    for j=1:nPerSideY
        x(end+1,:) = [i j];
    end
end

%plot(x(:,1),x(:,2),'.');
%hold  on
edges = zeros(size(x,1));
special_edges = zeros(size(x,1));
for i=1:length(edges)
    for j=1:length(edges)
        if sum(abs(x(i,:)-x(j,:)))==1
            edges(i,j)=1;
        end
    end
end

edges = edges-diag(diag(edges));
edges = double(triu(edges));
nVertices = size(edges,1);

if interact_sign~=4
  s_single = 2*(rand(nVertices,1)-0.5)*field_strength;
else
  s_single = gaus_mix(nVertices,[-1 0 1],[0.1 0.1 0.1],[1 0 1]/2);
  s_single = s_single*field_strength;
end

G.Fields = s_single;

switch interact_sign
 case 1 % Positive
  s_pair = (rand(nVertices).*edges);
 case 2
   s_pair = (-rand(nVertices).*edges);
 case 3
   s_pair = (2*(rand(nVertices)-0.5).*edges);
 case 4
   m = gaus_mix(nVertices^2,[-1 0 1],[0.1 0.1 0.1],[1 0 1]/2);
   s_pair = reshape(m,nVertices,nVertices).*edges;
end
s_pair = s_pair*interact_strength;

s_pair = s_pair+s_pair';
%s_pair(s_pair~=0) = exp(s_pair(s_pair~=0));
edges = edges+special_edges;

G.ArrWeights = s_pair;
edges = triu(edges)+triu(edges)';
G.Edges = edges;
G.List=[];
for i=1:size(s_pair,1)
    for j=1:size(s_pair,1)
        G.Weights(i,j)= [s_pair(i,j)];
        if s_pair(i,j) & i<j
            G.List(end+1,:) = [i j s_pair(i,j)];
        end
    end
end

G.FaceEdges={};
G.FaceVertices={};
List = G.List;

% Generate structures for Martin's JT
EdgeMap = zeros(size(edges));
for i=1:size(List,1)
  EdgeMap(List(i,1),List(i,2)) = i;
end
EdgeMap = EdgeMap+triu(EdgeMap)';
nVertices = size(x,1);
nEdges = size(List,1);
mu(1:nVertices) = G.Fields;
mu(nVertices+1:nVertices+nEdges) = List(:,3);
G.EdgeMap = EdgeMap;
G.mu = mu;