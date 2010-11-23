function qq = additem(qq,xx);
% qq = additem(qq,xx)
% adds data item xx into component qq.

qq.nn = qq.nn + 1;
qq.rr = qq.rr + 1;
qq.vv = qq.vv + 1;
qq.CC = cholupdate(qq.CC,xx);
qq.XX = qq.XX + xx;

