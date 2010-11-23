function pi = mean(qq);
% pi = mean(qq)
% returns the mean of pi given the data items in the component.

pi = (qq.aa + qq.mi) / (qq.aa*qq.dd + qq.mm);
