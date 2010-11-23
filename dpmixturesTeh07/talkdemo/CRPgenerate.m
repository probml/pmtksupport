% generates an assignment of customers to tables in CRP

aa = 10; % alpha
dd = 0; % discount
NN = 10000; % number of data points

nn = []; % row vector of number of customers per table
kk = zeros(1,NN); % k(i) = table number that customer i sat at 
KK = zeros(1,NN); % K(i) = number of occupied tables after first i customers

curK = 0;
for ii = 1:NN
  pp = [nn-dd aa+curK*dd]; % prop to prob of sitting at each table or new one
  kk(ii) = sum(rand(1)*sum(pp) > cumsum(pp))+1; % sample table number
  if kk(ii) <= length(nn)
    % occupied table
    nn(kk(ii)) = nn(kk(ii)) + 1;
  else
    % new table
    nn(kk(ii)) = 1;
    curK = curK + 1;
  end
  KK(ii) = curK;
end  
plot(kk,'.');
hold on
plot(KK,'-');
hold off
