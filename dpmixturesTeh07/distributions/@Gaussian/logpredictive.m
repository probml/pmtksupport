function ll = logpredictive(qq,xx);
% ll = logpredictive(qq,xx)
% log predictive probability of xx given other data items in the component
% log p(xx|x_1,...,x_n)

ll =   ZZ(qq.dd,qq.nn+1,qq.rr+1,qq.vv+1,cholupdate(qq.CC,xx),qq.XX+xx) ...
     - ZZ(qq.dd,qq.nn  ,qq.rr  ,qq.vv  ,           qq.CC    ,qq.XX   );
