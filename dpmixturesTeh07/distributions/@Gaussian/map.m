function [mu,sigma] = map(qq);
% [mu,sigma] = map(qq)
% Returns MAP estimate for mean and covariance of data items in component qq.

mu = qq.XX/qq.rr;
CC = cholupdate(qq.CC,qq.XX/sqrt(qq.rr),'-');
sigma = CC'*CC/(qq.vv-qq.dd-1);
