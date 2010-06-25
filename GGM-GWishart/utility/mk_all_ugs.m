function Gs = mk_all_ugs(N)
%
%Gs = mk_all_ugs(N)
%
% Generate all undirected graphs on N variables and output in 
% the cell array Gs.

nedges    = nchoosek(N,2);
m         = 2^nedges;
ind       = dec2bin(0:(m-1))-'0';  % all bit vectors
[foo,ord] = sort(sum(ind,2),'ascend');
ind       = ind(ord,:);
ut        = triu(ones(N),1)>0;

Gs = {};
for i=1:m
  G = zeros(N,N);
  G(ut) = ind(i,:);
  G = (G + G')>0;
  Gs{i} = G;
end
