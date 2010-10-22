% Version: Oct. 2004
% Author : Chen Yanover
% Overloading size function for the sparse_cell class 
function a = size(sc,d)

if nargin<2
    a = [sc.n, sc.m];
elseif d==1, a = sc.n;
elseif d==2, a = sc.m; 
else error('sparse_cell size: wrong usage'); end