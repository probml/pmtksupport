% Version: Oct. 2004
% Author : Chen Yanover
% Overloading display function for the sparse_cell class 
function display(sc)

disp(' ');
disp([inputname(1) ' = ']);
disp(' ');
disp(sprintf('     sparse_cell [%d, %d]', sc.n, sc.m));