% Version: Oct. 2004
% Author : Chen Yanover
% Overloading subassign function for the sparse_cell class 
function a = subsasgn(a,s,b)

switch s.type
case '{}',
    x = s.subs{1}; 
    if length([s.subs{:}])==1,
        y=1;
    else, y = s.subs{2}; end 
    if x<=a.n & y<=a.m,
        ind = a.indMat(x,y); 
    else 
        ind = 0; 
        a.n = max(a.n, x);
        a.m = max(a.m, y);
    end 
    
    if ind, a.cells{ind} = b;
    else 
        a.cells{end+1} = b; 
        a.indMat(x,y) = length(a.cells);
    end
case '()',
    nelem_ = length(b(:));
    a.cells(end+1:end+nelem_) = b(:)';
    a.indMat(s.subs{:}) = reshape(a.nelem+1:a.nelem+nelem_, size(b));
    [a.n, a.m] = size(a.indMat);
    a.nelem = a.nelem + nelem_;
end