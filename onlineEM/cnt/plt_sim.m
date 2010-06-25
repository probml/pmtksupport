% Script to plot the results of ex_sim
%
% $Id: plt_sim.m,v 1.1 2006/04/06 08:36:52 oKp Exp $

if (0)
  % Save the files produced by ex_sim
  R = input('save? ', 's');
  if (length(R) == 0)
    eval(['save ex_sim' int2str(NS) ' rate_rem wght_rem rate_titt ' ...
    'wght_titt']);
  end
end


% Plot
N = 20;
val = [10 50 100 500 1000 5000 10000 50000 100000];
WGHT_REM = zeros(N, length(val));
RATE1_REM = zeros(N, length(val));
RATE2_REM = zeros(N, length(val));
WGHT_TITT = zeros(N, length(val));
RATE1_TITT = zeros(N, length(val));
RATE2_TITT = zeros(N, length(val));
for i = 1:N
  load(['data/ex_sim' int2str(i)]);
  WGHT_REM(i,:) = wght_rem(val,1)';
  RATE1_REM(i,:) = rate_rem(val,1)';
  RATE2_REM(i,:) = rate_rem(val,2)';
  WGHT_TITT(i,:) = wght_titt(val,1)';
  RATE1_TITT(i,:) = rate_titt(val,1)';
  RATE2_TITT(i,:) = rate_titt(val,2)';
end
clf;
subplot(3,2,1);
boxplot(WGHT_REM);
ylabel('\pi_1');
xlabel('');
axis([0.5 9.5 0.5 1.1]);
subplot(3,2,2);
boxplot(WGHT_TITT);
ylabel('');
xlabel('');
axis([0.5 9.5 0.5 1.1]);

subplot(3,2,3);
boxplot(RATE1_REM);
ylabel('\lambda_1');
xlabel('');
axis([0.5 9.5 0.4 2.5]);
subplot(3,2,4);
boxplot(RATE1_TITT);
ylabel('');
xlabel('');
axis([0.5 9.5 0.4 2.5]);

subplot(3,2,5);
boxplot(RATE2_REM);
ylabel('\lambda_2');
xlabel('');
axis([0.5 9.5 0.2 7]);
subplot(3,2,6);
boxplot(RATE2_TITT);
ylabel('');
xlabel('');
axis([0.5 9.5 0.2 7]);

print -depsc plt_sim

save plt_sim val WGHT_REM RATE1_REM RATE2_REM WGHT_TITT RATE1_TITT RATE2_TITT;
