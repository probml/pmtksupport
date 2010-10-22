% Version: Oct. 2004
% Author : Chen Yanover
% Overloading subref function for the sparse_cell class 
function b = subsref(a,s)

s1 = s(1);
switch s1.type
case '()'
    im = subsref(a.indMat, s1);
    [b.n, b.m] = size(im);   
    nonEmpty = find(im);
    b.nelem = length(nonEmpty);
    c = a.cells(im(nonEmpty));
    im(nonEmpty) = 1:length(nonEmpty);
    b.indMat = im; 
    b.cells = c;
    b = class(b, 'sparse_cell');
case '{}'
    s1_ = s1; s1_.type = '()';
    im = subsref(a.indMat, s1_);
    if length(im(:))>1,
        error('Not yet supported'); end
    if im, b{1} = a.cells{im};
    else b{1} = []; end
    b = b{:};
otherwise
    error;
end

s = s(2:end);
if ~isempty(s),
    b = subsref(b,s); end
