
function [new_lambda,new_local] = negate_field(adj,lambda,local)

degrees = sum(adj,2);
nVars = size(adj,1);
new_lambda = cell(nVars);

% Divide the field

for i=1:nVars
  for j=1:nVars
    if adj(i,j)==0
      continue;
    end
    new_lambda{i,j} = -lambda{i,j};
  end
end

for i=1:nVars
  new_local{i} = -local{i}; 
end