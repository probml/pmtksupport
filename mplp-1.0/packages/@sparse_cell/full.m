% Version: Oct. 2004
% Author : Chen Yanover
% Overloading full function for the sparse_cell class 
function c = full(sc)

c = cell(sc.n, sc.m);
[x,y] = find(sc.indMat);

for i=[x(:),y(:)]',
    c(i(1),i(2)) = sc.cells(sc.indMat(i(1),i(2)));
end