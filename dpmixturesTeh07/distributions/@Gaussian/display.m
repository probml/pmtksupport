function display(qq);
% display GaussianGamma

if isequal(get(0,'FormatSpacing'),'compact')
  disp([inputname(1) ' = ']);
  disp(qq);
else
  disp(' ');
  disp([inputname(1) ' = ']);
  disp(' ');
  disp(qq);
end

