% minFunc
fprintf('Compiling minFunc files...\n');
mex -outdir minFunc minFunc/lbfgsC.c

% L1GeneralGroup
fprintf('Compiling L1GeneralGroup files...\n');
mex -outdir L1GeneralGroup/mex L1GeneralGroup/mex/auxGroupL2ProjectC.c
mex -outdir L1GeneralGroup/mex L1GeneralGroup/mex/auxGroupLinfProjectC.c
mex -outdir L1GeneralGroup/mex L1GeneralGroup/mex/projectRandom2C.c