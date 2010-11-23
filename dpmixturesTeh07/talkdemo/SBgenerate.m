% generates from stick-breaking construction
alpha = [5 5 5];
nn = 59;

for ii = 1:length(alpha)
subplot(length(alpha),1,ii);

beta = betarnd(1,alpha(ii), 1,nn);
neg = cumprod(1-beta);
pi = beta .* [1 neg(1:end-1)];

bar(1:nn,pi);
if ii == 1, title('stick-breaking weights pi'); end
if ii == length(alpha), xlabel(['stick indices']); end
ylabel(['\alpha = ' num2str(alpha(ii))]);
end
