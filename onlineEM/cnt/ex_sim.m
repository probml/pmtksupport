% Script to illustrate the behavior of the REM algorithm on simulated data
%
% $Id: ex_sim.m,v 1.3 2006/04/06 08:36:52 oKp Exp $

%%% % Fix random seeds
%%% %rand('state',7500833000)
%%% %randn('state',0142561451);

% Generate mixture data
T = 100000;
wght_sim = [0.8, 0.2];
rate_sim = [1, 2.5];
count = pm_gen(wght_sim, rate_sim, T);
% Starting values
wght_0 = wght_sim; % [0.5, 0.5];
rate_0 = rate_sim; % [1, 1.25];

% % Learn with EM
% disp('EM');
% Nit = 50;
% [wght_em, rate_em, logl_em] = pm_em(count, wght_0, rate_0, Nit);
% % Plots
% clf;
% plot(logl_em/T);

% Learn with REM
disp('REM');
gamma = (1:T).^(-0.5);
[wght_rem, rate_rem, logl_rem] = pm_rem(count, wght_0, rate_0, gamma);

% Learn with Titt
disp('Titt');
[wght_titt, rate_titt, logl_titt] = pm_titt(count, wght_0, rate_0, gamma);

% Plots
clf;
subplot(3,1,1)
plot(wght_rem(:,1));
hold on;
plot(wght_titt(:,1),'k');
%axis([0 T 0.5 1]);
ylabel('\pi_1');
plot(average(wght_rem(:,1),0.999,'exp'), 'r');
plot(average(wght_rem(:,1),1000,'mean'), 'm');

subplot(3,1,2)
plot(rate_rem(:,1));
hold on;
plot(rate_titt(:,1),'k');
%axis([0 T 0.5 1.5]);
ylabel('\lambda_1');
plot(average(rate_rem(:,1),0.999,'exp'), 'r');
plot(average(rate_rem(:,1),1000,'mean'), 'm');

subplot(3,1,3)
plot(rate_rem(:,2));
hold on;
plot(rate_titt(:,2),'k');
%axis([0 T 4.5 5.5])
ylabel('\lambda_2');
xlabel('Number of iterations');
plot(average(rate_rem(:,2),0.999,'exp'), 'r');
plot(average(rate_rem(:,2),1000,'mean'), 'm');
