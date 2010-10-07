clear all
close all

bugsFolder = 'C:\kmurphy\Programs\WinBUGS14';
%demoFolder = 'C:\kmurphy\Programs\WinBUGS14\Manuals\Tutorial';
demoFolder = 'C:\kmurphy\GoogleCode\matbugs\demos\seedsSimpleDemo';

% read data
% r n x1 x2
data = importdata(fullfile(demoFolder, 'seeds_data.txt'), '\t', 1);
data.data(:,1) = []; % strip off empty column
dataStruct = struct( ...
  'r', data.data(:,1), ...
  'n', data.data(:,2), ...
  'x1', data.data(:,3), ...
  'x2', data.data(:,4), ...
  'N', size(data.data,1));
 

  
Nchains = 2;

% We can use a gamma prior on precision tau
% or a uniform prior on sd sigma
useGamma = true;

% specify initial values by cutting and pasting from init file
S.alpha0 = 0; S.alpha1 = 0; S.alpha2 = 0; S.alpha12 = 0;
if useGamma
  S.tau = 1;
else
  S.sigma = 1;
end
S.b = zeros(1,21); 
initStructs(1) = S;

% Actually we had to modify the init file
% to make sigma be smaller
% We also reduce the size of alpha
% Otherwise bugs crashes!
S.alpha0 = 1; S.alpha1 = 1; S.alpha2 = 1; S.alpha12 = 1;
if useGamma
  S.tau = 0.1;
else
  S.sigma= 2;
end
S.b = [0.1, -0.2, 0.25, 0.11, -0.21, 0.3, -0.25, 0.15, -0.31, -0.1,...
0.1, 0.12, 0.2, -0.2, 0.4, -0.24, 0.14, 0.3, -0.2, 0.1, 0.05];
initStructs(2) = S;


if useGamma
  model = 'seeds_model.txt';
else
  model = 'seeds_unif_model.txt';
end
% run sampler
[samples, stats] = matbugs(dataStruct, ...
  fullfile(demoFolder, model), ...
  'init', initStructs, ...
  'nchains', 2, ...
  'view', 1, 'nburnin', 1000, 'nsamples', 2000, ...
  'thin', 1, ...
  'monitorParams', {'alpha0','alpha1','alpha2','alpha12'}, ...
  'Bugdir', bugsFolder);

