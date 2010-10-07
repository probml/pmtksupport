
         
% The data:

r = [10, 23, 23, 26, 17, 5, 53, 55, 32, 46, 10,...
    8, 10,   8, 23, 0,  3, 22, 15, 32, 3];
n = [39, 62, 81, 51, 39, 6, 74, 72, 51, 79, 13, 16,...
    30, 28, 45, 4, 12, 41, 30, 51, 7];
x1 = [0,   0,  0,   0,   0, 0,   0,   0,  0,   0,   0,...
    1,   1,   1,   1, 1,   1,  1,   1,   1, 1]; % seed type
x2 = [0,   0,  0,   0,   0, 1,   1,   1,  1,   1,   1,  0,...
    0,   0,   0, 0,   1,  1,   1,   1, 1]; % root extract
N = 21;

% some common input parameters you might with to change for all calls to
% winbugs:
burn = 1000; % number of iterations for burn in
samp = 2000; % number of iterations to collect data
bugsFolder = 'C:\kmurphy\Programs\WinBUGS14';
%bugsFolder = 'C:/Program Files/WinBUGS14';

% run the model using both seed and root extract type as factors
% as in the example in Winbugs
dataseed = struct('r',r,'n',n,'x1',x1,'x2',x2,'N',N);
init1 = struct('alpha0',0,'alpha1',0,'alpha2',0,'alpha12',0,'sigma',1);

[fullsamples, fullstats] = matbugs(dataseed, ...
		fullfile(pwd, 'seeds_unif_model.txt'), ...
		'init', init1, ...
		'nChains', 1, ...
		'view', 0, 'nburnin', burn, 'nsamples', samp, ...
		'thin', 1, 'DICstatus', 1, ...
        'monitorParams', {'alpha0','alpha1','alpha2','alpha12','sigma'}, ...
		'Bugdir', bugsFolder);
         
         
         
         
         
         
         
         
         
         
         
         
         
         