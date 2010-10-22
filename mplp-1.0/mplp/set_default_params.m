function params=set_default_params(params,fields,values,n)

for i=1:length(fields)
  if ~isfield(params,fields{i})
     eval(sprintf('params.%s= values(i);',fields{i}));	
  end
end

if params.file_type == 1
  if xor(isfield(params,'num_rows'), isfield(params,'num_cols'))
    error('One of params.num_rows or params.num_cols is missing.');
  elseif ~isfield(params,'num_rows')
    
    if sqrt(n) ~= round(sqrt(n))  
      error(['Number of variables is not a perfect square; Please set ' ...
             'params.num_rows or params.num_cols.']);
    end
  
    params.num_rows = sqrt(n);
    params.num_cols = sqrt(n);
  end
end