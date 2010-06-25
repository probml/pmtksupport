% Script to run EM on the toy example taken from flexmix (1 and 6 used to produce figures for
% the paper)
%
% $Id: ex_toy.m,v 1.6 2007/10/01 15:40:27 cappe Exp $

% 1. Simulate data
% Dimensions
n = 3000;
m = 2;
d = 3;
% Simulate data
U = 10*rand(1,n); % Easy one; U = 5+randn(1,n); % Hard one
W = 1 + round(rand(1,n));
% Regressions parameters
beta_t = [[0; 5; 0], [15; 10; -10]];
sigma2_t = [81 81];
w_t = [0.5 0.5];
Y = zeros(1,n);
Y(W == 1) = beta_t(2,1)*U(W == 1);
Y(W == 2) = beta_t(1,2)*ones(1, sum(W == 2)) + ...
  beta_t(2,2)*U(W == 2) + beta_t(3,2)*(U(W == 2)).^2/10;
Y = Y + sqrt(sigma2_t(1))*randn(1,n);
Z = [ones(1,n); U; U.^2/10];
% Plot
figure(1);
clf;
plot(U(W == 1), Y(W == 1), 'ok');
hold on;
plot(U(W == 2), Y(W == 2), 'xr');
xlabel('U');
ylabel('Z');
% graph('batch_data500');

% 2. Estimate information matrix (only for betas)
inf1 = zeros(d,d,n);
inf2 = zeros(d,d,n);
for i = 1:n
   [inf1_one, inf2_one] = reg_info(Y(i), Z(:,i), w_t, beta_t, sigma2_t);
   inf1(:,:,i) = inf1_one(:,:,1);
   inf2(:,:,i) = inf2_one(:,:,1);
end
% %%% n = 1000000
% >> I = inv(mean(inf1,3))
% I =
%
%    1.0e+03 *
%
%     2.2859   -0.9183    0.7588
%    -0.9183    0.4882   -0.4512
%     0.7588   -0.4512    0.4439
%    
% >> s = sqrt(diag(I))
%
% s =
%
%    47.8108
%    22.0954
%    21.0682
%
% >> I ./ (s*s')
%
% ans =
%
%     1.0000   -0.8692    0.7533
%    -0.8692    1.0000   -0.9693
%     0.7533   -0.9693    1.0000


% 3. Some starting values
w_0 = [0.5 0.5];
beta_0 = [[0; 4.9; 0], [30; 0; 0]];
sigma2_0 = [230, 230];

% 4. EM
nit = 100;
[w, beta, sigma2, logl] = reg_batch(Y, Z, w_0, beta_0, sigma2_0, nit);
% Plot
figure(2);
clf;
plot(logl);
hold on;
% For comparison purposes (likelihood at true param.)
[w_t, beta_t, sigma2_t, logl_t] = ...
    reg_batch(Y, Z, w_t, beta_t, sigma2_t, 1);
plot([0, nit+1], logl_t(1)*ones(1,2), 'r');

% 5. On-line EM
gam = (1:n).^(-0.6);
[w, beta, sigma2, logl] = reg_online(Y, Z, w_0, beta_0, sigma2_0, gam);
clf;
plot(squeeze(beta(1,1,:)),'k');
hold on;
plot(squeeze(beta(2,1,:)),'k:');
plot(squeeze(beta(3,1,:)),'k-.');


% 6. Plot some trajectories
gam = (1:n).^(-1);
[w, beta, sigma2, logl] = reg_online(Y, Z, w_0, beta_0, sigma2_0, gam);
clf;
for i = 1:3
  subplot(3,1,i);
  plot(squeeze(beta(i,2,:)),'--r');
  hold on;
  if (i == 1)
    axis([0 n -5 35]);
  elseif (i ==2)
    axis([0 n 0 20]);
  else
    axis([0 n -20 0]);
  end    
  ylabel(['\beta_2(' int2str(i) ')']);
end
gam = (1:n).^(-0.6);
[w, beta, sigma2, logl] = reg_online(Y, Z, w_0, beta_0, sigma2_0, gam);
for i = 1:3
  subplot(3,1,i);
  plot(squeeze(beta(i,2,:)),':k');
end
% With averaging
nav = 1000;
for i = 1:3
  subplot(3,1,i);
  t = squeeze(beta(i,2,:));
  t(nav+1:n) = cumsum(t(nav+1:n))./(1:n-nav)';
  plot(t,'-b');
end
xlabel('samples');
%graph('batch_traj_3000a1000')
