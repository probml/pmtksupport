function display(qq);
% display Multinomial

if isequal(get(0,'FormatSpacing'),'compact')
  disp([inputname(1) ' = ']);
  disp(qq);
else
  disp(' ');
  disp([inputname(1) ' = ']);
  disp(' ');
  disp(qq);
end

