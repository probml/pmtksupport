function [d] = lbfgsC(g, s, y, Hdiag)
%% Call lbfgs if lbfgsC.mex* is not available. 
d = lbfgs(g, s, y, Hdiag); 
end