% Code to demonstrate DP mixture models, and also used in the practical course.
% Written by Yee Whye Teh.
%
% Subdirectories:
% distributions         - implements Gaussian and Multinomial distributions
% talkdemo              - The demos used in the talk
% dpmixture             - DP mixture models
% finitemixture         - finite mixture models
% utilities             - some utilities
% nipsdata              - NIPS data
%
% help distributions or dpmixture or finitemixture for more information
% 
% run initpath in matlab/ directory to initialize paths in matlab
%
% Demos:
% talkdemo/GPgenerate           - generates from Gaussian process
% talkdemo/DPgenerate           - generates from Dirichlet process
% talkdemo/DPposterior          - generates from posterior Dirichelt process
% talkdemo/SBgenerate           - generates from stick-breaking construction
% talkdemo/DPmixgenerate        - generates from DP mixture model
% talkdemo/CRPgenerate          - generates from Chinese restaurant process
% finitemixture/fm_demo1d       - inference in finite mixture of Gaussians in 1D
% finitemixture/fm_demo2d       - inference in finite mixture of Gaussians in 2D
% dpmixture/dpm_demo1d          - inference in DP mixture of Gaussians in 1D
% dpmixture/dpm_demo2d          - inference in DP mixture of Gaussians in 2D
% dpmixture/dpm_demonips        - DP mixture for clustering NIPS papers
