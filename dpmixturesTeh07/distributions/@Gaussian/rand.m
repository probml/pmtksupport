function [mu,sigma] = rand(qq);
% [mu, sigma] = rand(qq)
% generates random mean mu and covariance sigma from the posterior of
% component parameters given data items in component.

CC = cholupdate(qq.CC,qq.XX/sqrt(qq.rr),'-')\eye(qq.dd);

sigma = iwishrnd(CC,qq.vv,CC);
mu = mvnrnd(qq.XX'/qq.rr, sigma/qq.rr)';
