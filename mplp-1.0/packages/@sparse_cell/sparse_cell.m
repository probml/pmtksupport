% Version: Oct. 2004
% Author : Chen Yanover
% 
% sparse_cell: creates a sparse cell
%     S = sparse_cell(a) converts a sparse or full cell to sparse form by
%     squeezing out any empty elements.
%  
%     S = sparse_cell(i,j,s,m,n) uses the rows of [i(:),j(:),s(:)] to generate
%     an m-by-n sparse cell.  The two integer index vectors, i and j, and the 
%     entries cell array s, all have the same length. 
%  
%     There are several simplifications of this five argument call:
%     1. S = sparse_cell(i,j,s) uses m = max(i) and n = max(j).
%     2. S = sparse_cell(m,n) abbreviates sparse_cell([],[],[],m,n).  This
%        generates the ultimate sparse cell, an m-by-n all empty cell array.
%     3. S = sparse_cell creates an empty sparse_cell 
function sc = sparse_cell(i,j,s,n,m)

% class members: 
%   n,m     - cell matrix dimensions
%   nelem   - number of non empty elements
%   indMat  - indices into the cell array
%   cells   - cell array
switch nargin,
case 0,     % creates an empty sparse cell
    sc.n = 0;
    sc.m = 0;
    sc.nelem = 0;
    sc.indMat = sparse([]);
    sc.cells = {};
    sc = class(sc, 'sparse_cell');
case 1,     % converts a sparse or full cell to sparse form
    if isa(i, 'sparse_cell') sc=i;
    elseif iscell(i),
        [sc.n, sc.m] = size(i);
        nelem = 0;
        indMat = sparse(sc.n, sc.m);
        cells = {};
        for n_=1:sc.n,
            for m_=1:sc.m,
                if ~isempty(i{n_,m_}),
                    nelem = nelem+1;
                    indMat(n_,m_) = nelem;
                    cells{end+1} = i{n_,m_};
                end
            end
        end
        sc.nelem = nelem;
        sc.indMat = indMat;
        sc.cells = cells;
        sc = class(sc, 'sparse_cell');
    else, error('sparse_cell: wrong usage'); 
    end
case {2,3,5},
    if nargin==2, 
        n=i; m=j;
        i=[]; j=[]; s={}; 
    elseif nargin==3,
        n = max(i); 
        m = max(j);
    elseif nargin==4, error('sparse_cell: wrong usage'); 
    end
    
    if ~iscell(s), 
        error('sparse_cell: conversion to cell from double is not possible'); end
    
    i=i(:); j=j(:); s=s(:); 
    if length(i)~=length(j) | length(j)~=length(s),
        error('sparse_cell: vectors must be the same lengths'); end
    
    sc.n = n; 
    sc.m = m;
    sc.nelem = length(i);
    sc.indMat = sparse(i, j, 1:sc.nelem, sc.n, sc.m);
    sc.cells = s;
    sc = class(sc, 'sparse_cell');
otherwise, error('sparse_cell: wrong usage');
end