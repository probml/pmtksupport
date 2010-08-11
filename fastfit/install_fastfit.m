cd(fileparts(which(mfilename))); 
mex pochhammer.c ../lightspeed2.3/util.obj -I../lightspeed2.3
mex di_pochhammer.c ../lightspeed2.3/util.obj -I../lightspeed2.3
%mex tri_pochhammer.c ../lightspeed/util.obj -I../lightspeed
%mex s_derivatives.c ../lightspeed/util.obj -I../lightspeed
