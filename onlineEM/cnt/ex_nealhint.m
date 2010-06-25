% Script to illustrate the behavior of the REM algorithm on simulated data
%
% $Id: ex_nealhint.m,v 1.1 2006/07/12 16:51:31 cappe Exp $

%%% % Fix random seeds
rand('state',7500833000)
randn('state',0142561451);

% Generate mixture data
T = 1000;
wght_sim = [0.8, 0.2];
rate_sim = [1, 3];
count = pm_gen(wght_sim, rate_sim, T);
% Starting values
wght_0 = [0.5, 0.5];
rate_0 = [1, 1.5];

% Learn with EM
disp('EM');
Nit = 150;
[wght_em, rate_em, logl_em] = pm_em(count, wght_0, rate_0, Nit);
% Plots
clf;
plot(logl_em);

% Learn with Neal-Hinton
disp('Neal-Hinton');
[wght_nh, rate_nh, logl_nh] = pm_nealhint(count, wght_0, rate_0, Nit);
[wght_nh5, rate_nh5, logl_nh5] = pm_nealhint(count, wght_0, rate_0, Nit, 5);
% Plot
hold on;
plot(logl_nh, 'r--');
plot(logl_nh5, 'r-.');

% Lean with REM
disp('REM');
gamma = (1:Nit*T).^(-0.5);
count_l = reshape((ones(Nit,1)*count)', 1, T*Nit);
[wght_rem, rate_rem, logl_rem] = pm_rem(count_l, wght_0, rate_0, gamma);
logl_rem = zeros(1, Nit);
for i = 1:Nit
  [wght_tmp, rate_tmp, logl_tmp] = pm_em(count, ...
    wght_rem(1+(i-1)*T, :), rate_rem(1+(i-1)*T, :), 1);
  logl_rem(i) = logl_tmp;
end
plot(logl_rem, 'k:');

legend('EM', 'NH', 'NH*5', 'REM', 4);
