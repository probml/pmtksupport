function [tt,cc] = autocorr(ss);

% computes the autocorrelation time of a sequence ss.

ss = ss - mean(ss);
cc = xcorr(ss,'unbiased');
cc = cc(end/2+1.5:end)/cc(end/2+.5);
mm = min(find(cc<0))-1;
tt = 1+2*sum(cc(1:mm));
