function disp(qq);
% disp Multinomial

[ii jj mi] = find(qq.mi);
disp(['Multinomial: d=' num2str(qq.dd) ...
                  ' a=' num2str(qq.aa*qq.dd) ...
                  ' n=' num2str(qq.nn) ...
                  ' m=' num2str(qq.mm)]);
if qq.mm>0
  disp(sprintf('%d:%d ',cat(1,ii',mi')));
end
