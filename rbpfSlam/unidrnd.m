function r = unidrnd(n,mm,nn)


if nargin == 1
    [errorcode rows columns] = rndcheck(1,1,n);
elseif nargin == 2
    [errorcode rows columns] = rndcheck(2,1,n,mm);
elseif nargin == 3
    [errorcode rows columns] = rndcheck(3,1,n,mm,nn);
else
    error('Requires at least one input argument.'); 
end

if errorcode > 0
    error('Size information is inconsistent.');
end

r = ceil(n .* rand(rows,columns));

% Fill in elements corresponding to illegal parameter values
if prod(size(n)) > 1
    r(n < 0 | round(n) ~= n) = NaN;
elseif n < 0 | round(n) ~= n
    r(:) = NaN;
end
