function N = pm_chk(wght, rate)

%pm_chk	  Checks the parameters of a Poisson mixture and returns its dimension.
%         Use: N = pm_chk(wght,rate).

% H2M/cnt Toolbox, Version 2.0
% Olivier Cappé, 30/12/97 - 30/12/97
% ENST Dpt. TSI / LTCI (CNRS URA 820), Paris

% Constants
% Max deviation from 1 (beware parameters may have been computed in float)
MAX_DEV = 1e-6;

% Check input arguments
error(nargchk(2, 2, nargin));
% Number of mixture components
N = length(wght);
if (any(wght < 0) | any(wght > 1))
  error('Mixture weights must be between 0 and 1.');
end
if (abs(sum(wght)-1) > MAX_DEV)
  error('Mixture weigths are not normalized.');
end
if (length(rate) ~= N)
  error('Vectors of mixture weigths and intensity must have the same size.');
end
if (any(rate) <= 0)
  error('Poisson rates must be positive.');
end
