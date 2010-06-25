function aest = average(est, param, mode);

%average
%         Use: average(est, param[, mode]);
%
% $Id: average.m,v 1.1 2006/04/04 17:05:39 oKp Exp $

if (nargin < 3)
  mode = 'mean';
end

aest = zeros(size(est));
l = length(est);
est = est(:);

if (strcmp(mode,'mean'))
  last = param;
  aest(1:last) = est(1:last);
  aest(last+1:l) = cumsum(est(last+1:l))./(1:(l-last))';
elseif (strcmp(mode,'exp'))
  % Exponential smoothing
  lambda = param;
  if ((lambda < 0) | (lambda > 1))
    error('Smoothing parameter must be between 0 and 1');
  end
  aest(1) = est(1);
  for i = 2:length(est);
    aest(i) = (1-lambda)*est(i) + lambda*aest(i-1);
  end
end
