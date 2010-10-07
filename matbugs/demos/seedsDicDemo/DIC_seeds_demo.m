% DIC demo
% Andrew Jackson 10 April  2006
% Social and Public Health Sciences Unit, Medical Research Council
% Glasgow G12 8RZ, UK.
% andrew@msoc.mrc.gla.ac.uk

% Uses the seeds example from the WinBUGS manual...
% The text from this example file is pasted below:
% This mfile runs 4 versions of this example:
% 1) with only seed type as an explanatory factor
% 2) with only root extract type as an explanatory factor
% 3) with both as factors but without an interation between them
% 4) with both as factors with an interation term
% Model fit is assessed through the DIC values. 
%
% TAKEN FROM WINBUGS MANUAL
% Spiegelhalter, D. J., Thomas, A., Best, N. G. & Lunn, D. 
%   WinBUGS Version 1.4.1 User Manual 
%   (Medical Research Council Biostatistics Unit, Cambridge, 2004).
% 
% % This example is taken from Table 3 of Crowder (1978), and concerns the 
% proportion of seeds that germinated on each of 21 plates arranged 
% according to a 2 by 2 factorial layout by seed and type of root extract. 
% The data are shown below, where ri and ni are the number of germinated 
% and the total number of seeds on the i th plate, i =1,...,N. 
% These data are also analysed by, for example, Breslow and Clayton (1993). 
% 
% 		seed O. aegyptiaco 75			seed O. aegyptiaco 73
% 	Bean            Cucumber        Bean            Cucumber
% 	r	n	r/n     r	n	r/n     r   n	r/n     r	n	r/n
% 	_________________________________________________________________
% 	10	39	0.26	5	6	0.83	8	16	0.50	3	12	0.25
% 	23	62	0.37	53	74	0.72	10	30	0.33	22	41	0.54	
% 	23	81	0.28	55	72	0.76	8	28	0.29	15	30	0.50
% 	26	51	0.51	32	51	0.63	23	45	0.51	32	51	0.63
% 	17	39	0.44	46	79	0.58	0	4	0.00	3	7	0.43
%                   10	13	0.77	
% 
% The model is essentially a random effects logistic, allowing for 
% over-dispersion.  If pi is the probability of germination on 
% the i th plate,  we assume
% 
% 	ri  ~  Binomial(pi, ni)
% 	
% 	logit(pi) = a0 + a1x1i + a2x2i + a12x1ix2i + bi
% 	
% 	bi  ~ Normal(0, t)
% 	
% where x1i  , x2i   are the seed type and root extract of the i th plate, 
% and an interaction term a12x1ix2i   is included.   a0 , a1 ,  a2 ,   a12  
% are given independent "noninformative" priors.

         
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

%{
% run the model using only the seed type as a factor
dataseed = struct('r',r,'n',n,'x1',x1,'N',N);
init1 = struct('alpha0',0,'alpha1',0,'sigma',1);


[x1samples, x1stats] = matbugs(dataseed, ...
		fullfile(pwd, 'seeds_x1.txt'), ...
		'init', init1, ...
		'nChains', 1, ...
		'view', 0, 'nburnin', burn, 'nsamples', samp, ...
		'thin', 1, 'DICstatus', 1, ...
        'monitorParams', {'alpha0', 'alpha1', 'sigma'}, ...
		'Bugdir', bugsFolder);
    
        
% run the model using only the root extract type as a factor
dataseed = struct('r',r,'n',n,'x2',x2,'N',N);
init1 = struct('alpha0',0,'alpha2',0,'sigma',1);

[x2samples, x2stats] = matbugs(dataseed, ...
		fullfile(pwd, 'seeds_x2.txt'), ...
		'init', init1, ...
		'nChains', 1, ...
		'view', 0, 'nburnin', burn, 'nsamples', samp, ...
		'thin', 1, 'DICstatus', 1, ...
        'monitorParams', {'alpha0','alpha2','sigma'}, ...
		'Bugdir', bugsFolder);
    
% run the model using both seed and root extract type as factors
% but with no covariate term.
dataseed = struct('r',r,'n',n,'x1',x1,'x2',x2,'N',N);
init1 = struct('alpha0',0,'alpha1',0,'alpha2',0,'sigma',1);

[nocovarsamples, nocovarstats] = matbugs(dataseed, ...
		fullfile(pwd, 'seeds_nocovar.txt'), ...
		'init', init1, ...
		'nChains', 1, ...
		'view', 0, 'nburnin', burn, 'nsamples', samp, ...
		'thin', 1, 'DICstatus', 1, ...
        'monitorParams', {'alpha0','alpha1','alpha2','sigma'}, ...
		'Bugdir', bugsFolder);
    
%}
% run the model using both seed and root extract type as factors
% as in the example in Winbugs
dataseed = struct('r',r,'n',n,'x1',x1,'x2',x2,'N',N);
init1 = struct('alpha0',0,'alpha1',0,'alpha2',0,'alpha12',0,'sigma',1);

[fullsamples, fullstats] = matbugs(dataseed, ...
		fullfile(pwd, 'seeds_full_model.txt'), ...
		'init', init1, ...
		'nChains', 1, ...
		'view', 0, 'nburnin', burn, 'nsamples', samp, ...
		'thin', 1, 'DICstatus', 1, ...
        'monitorParams', {'alpha0','alpha1','alpha2','alpha12','sigma'}, ...
		'Bugdir', bugsFolder);
         
disp('DIC for seed type only')         
%disp(x1stats.DIC.total)
disp(x1stats.DIC.S.DIC)
disp('DIC for root extract type only')
disp(x2stats.DIC.S.DIC)
disp('DIC allowing for both seed and root type WITHOUT an interaction')
disp(nocovarstats.DIC.S.DIC)
disp('DIC allowing for both seed and root type with an interaction')
disp(fullstats.DIC.S.DIC)
sprintf('... including both explanatory factors with an interaction term\nincreases model fit as determined by lower DIC values')

         
         
         
         
         
         
         
         
         
         
         
         
         
         
         
         
         
         