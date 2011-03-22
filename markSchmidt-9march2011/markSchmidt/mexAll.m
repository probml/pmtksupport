% minFunc
fprintf('Compiling minFunc files...\n');
mex -outdir minFunc minFunc/mcholC.c
mex -outdir minFunc minFunc/lbfgsC.c

% KPM
fprintf('Compiling KPM files...\n');
mex -IKPM -outdir KPM KPM/max_mult.c

% DAGlearn
fprintf('Compiling DAGlearn files...\n');
mex -IDAGlearn/ancestorMatrix -outdir DAGlearn/ancestorMatrix DAGlearn/ancestorMatrix/ancestorMatrixAddC_InPlace.c
mex -IDAGlearn/ancestorMatrix -outdir DAGlearn/ancestorMatrix DAGlearn/ancestorMatrix/ancestorMatrixBuildC.c

% GroupL1
fprintf('Compiling L1GeneralGroup files...\n');
mex -outdir L1GeneralGroup/mex L1GeneralGroup/mex/projectRandom2C.c
mex -outdir L1GeneralGroup/mex L1GeneralGroup/mex/auxGroupLinfProjectC.c
mex -outdir L1GeneralGroup/mex L1GeneralGroup/mex/auxGroupL2ProjectC.c

% UGM
fprintf('Compiling UGM files...\n');
mex -IUGM/mex -outdir UGM/mex UGM/mex/UGM_makeNodePotentialsC.c
mex -IUGM/mex -outdir UGM/mex UGM/mex/UGM_makeEdgePotentialsC.c
mex -IUGM/mex -outdir UGM/mex UGM/mex/UGM_MRFLoss_subC.c
mex -IUGM/mex -outdir UGM/mex UGM/mex/UGM_updateGradientC.c
mex -IUGM/mex -outdir UGM/mex UGM/mex/UGM_PseudoLossC.c
mex -IUGM/mex -outdir UGM/mex UGM/mex/UGM_Decode_ICMC.c
mex -IUGM/mex -outdir UGM/mex UGM/mex/UGM_Loss_subC.c
mex -IUGM/mex -outdir UGM/mex UGM/mex/UGM_Infer_LBPC.c
mex -IUGM/mex -outdir UGM/mex UGM/mex/UGM_Infer_ExactC.c
mex -IUGM/mex -outdir UGM/mex UGM/mex/UGM_Sample_GibbsC.c
mex -IUGM/mex -outdir UGM/mex UGM/mex/UGM_Decode_ExactC.c
mex -IUGM/mex -outdir UGM/mex UGM/mex/UGM_Infer_ChainC.c
mex -IUGM/mex -outdir UGM/mex UGM/mex/UGM_Infer_MFC.c
mex -IUGM/mex -outdir UGM/mex UGM/mex/UGM_Infer_TRBPC.c

% Newer UGM stuff
mex -IUGM/mex -outdir UGM/mex UGM/mex/UGM_CRF_makePotentialsC.c
mex -IUGM/mex -outdir UGM/mex UGM/mex/UGM_LogConfigurationPotentialC.c
mex -IUGM/mex -outdir UGM/mex UGM/mex/UGM_CRF_NLLC.c

% Ewout
fprintf('Compiling Ewout files...\n');
mex -IForeign/Ewout -outdir Foreign/Ewout Foreign/Ewout/projectBlockL1.c Foreign/Ewout/oneProjectorCore.c Foreign/Ewout/heap.c
mex -IForeign/Ewout -outdir Foreign/Ewout Foreign/Ewout/projectBlockL2.c

% crfChain
fprintf('Compiling crfChain files...\n');
mex -outdir crfChain/mex crfChain/mex/crfChain_makePotentialsC.c
mex -outdir crfChain/mex crfChain/mex/crfChain_inferC.c
mex -outdir crfChain/mex crfChain/mex/crfChain_lossC2.c