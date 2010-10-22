% Version: Oct. 2004
% Author : Chen Yanover
% Overloading iscell function for the sparse_cell class 
function t = iscell(a)

t = isa(a, 'sparse_cell');
