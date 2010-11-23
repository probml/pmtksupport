function disp(qq);
% disp GaussianWishart

[mu sigma] = map(qq);
disp(['GaussianWishart: d=' num2str(qq.dd) ...
                      ' n=' num2str(qq.nn) ...
                      ' r=' num2str(qq.rr) ...
                      ' v=' num2str(qq.vv)]);
disp(' mu=');
disp(mu);
disp(' sigma=');
disp(sigma);
