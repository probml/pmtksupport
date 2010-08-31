function [d] = lbfgsC(g,s,y,Hdiag)
%% See lbfgs
% This is called if a lbfgsC.mex* file cannot be found your system. 
d = lbfgs(g, s, y, Hdiag); 

end