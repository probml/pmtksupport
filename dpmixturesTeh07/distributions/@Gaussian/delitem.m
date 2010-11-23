function qq = delitem(qq,xx);
% qq = delitem(qq,xx)
% deletes data item xx from component qq.

qq.nn = qq.nn - 1;
qq.rr = qq.rr - 1;
qq.vv = qq.vv - 1;
qq.CC = cholupdate(qq.CC,xx,'-');
qq.XX = qq.XX - xx;

