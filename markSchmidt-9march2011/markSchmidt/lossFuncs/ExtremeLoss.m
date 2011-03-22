function [nll,g,H] = ExtremeLoss(w,X,y)
% w(feature,1)
% X(instance,feature)
% y(instance,1)
%
% Binary response with complementary log-log link

Xw = X*w;
p_1 = 1 - exp(-exp(Xw));
nll = -sum(log(p_1(y==1)))-sum(log(1-p_1(y==-1)));

tmp = exp(Xw - exp(Xw));
g = -X(y==1,:)'*(tmp(y==1)./p_1(y==1)) - X(y==-1,:)'*(-tmp(y==-1)./(1-p_1(y==-1)));

