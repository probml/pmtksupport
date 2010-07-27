function [count, label] = pm_gen(wght, rate, T)

%pm_gen	  Simulates data from a Poisson mixture.
%         Use: [count,label] = pm_gen(wght,rate,T) where count is an
%         array of length T which contains the simulated data and
%         label contains the corresponding simulated mixture indicators.
%         Requires Matlab's Statistics toolbox or GNU Octave.

% H2M/cnt Toolbox, Version 2.0
% Olivier Cappé, 30/12/97 - 22/08/2001
% ENST Dpt. TSI / LTCI (CNRS URA 820), Paris

% Needed functions
if (~((exist('randindx') == 2) | (exist('randindx') == 3)))
  error('Function randindx (from h2m toolbox) is missing.');
end
% Check that we have something to simulate from Poisson
if (exist('poisson_rnd') == 2)
  S_VER = 0; % Octave
elseif (exist('poissrnd') == 2)
  S_VER = 1; % Matlab + statistics toolbox
else
  error(['Either poissrnd (Matlab''s Statistics toolbox) or\n' ...
         'poisson_rnd (Octave >= 2.014) are required (but read\n' ...
	 'the H2M doc to see how to get around the problem.)']);
end

% Inputs arguments
error(nargchk(3, 3, nargin));
N = pm_chk(wght, rate);

% First simulate labels
label = randindx(wght, T);

% Then the Poisson data
if (S_VER == 0)
  count = poisson_rnd(rate(label));
else
  count = poissrnd(rate(label));
end
