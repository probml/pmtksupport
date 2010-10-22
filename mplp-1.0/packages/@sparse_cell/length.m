% Version: Oct. 2004
% Author : Chen Yanover
% Overloading length function for the sparse_cell class 
function k = length(sc)

k = max(sc.m, sc.n);