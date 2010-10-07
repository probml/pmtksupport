% An example illustrating how to call matbugs.m
% We use the schools data/model from p140/ p592 of Gelman's book
% Written by Maryam Mahdaviani, August 2005
% Modified by Kevin Murphy (murphyk@cs.ubc.ca), 11 October 2007

dataStruct = struct('J', 8, ...
                    'y', [28, 8, -3, 7, -1, 1, 18, 12], ...
                   'sigma_y', [15, 10, 16, 11, 9, 11, 10, 18]);

Nchains = 3;

% we initialize the params to the observed data values, but with decreasing
% confidence, as suggested on p593 of Gelman
clear initStructs
for i=1:Nchains
  S.theta = dataStruct.y;
  S.mu_theta = 0;
  S.sigma_theta = 10^i; % each chain becomes more over-dispersed
  initStructs(i) = S;
end

bugsFolder = 'C:\kmurphy\Programs\WinBUGS14';
[samples, stats] = matbugs(dataStruct, ...
		fullfile(pwd, 'schools_model.txt'), ...
		'init', initStructs, ...
		'view', 0, 'nburnin', 1000, 'nsamples', 500, ...
		'thin', 10, ...
		'monitorParams', {'theta', 'mu_theta', 'sigma_theta'}, ...
		'Bugdir', bugsFolder);

fprintf('Rhat\n');
stats.Rhat

fprintf('means\n');
stats.mean

% Traceplots
figure;
colors = 'rgb';
for j=1:8
  subplot(3,3,j); hold on
  for c=1:Nchains
    %plot(s(c).theta(:,j), colors(c));
    plot(samples.theta(c,:,j), colors(c));
  end
  title(sprintf('theta %d', j));
end
subplot(3,3,9); hold on
for c=1:Nchains
  plot(samples.mu_theta(c,:), colors(c));
end
title(sprintf('mu.theta'))

% Posterior summaries - kernel density estimation
figure;
for j=1:8
  subplot(3,3,j); hold on
  %dat = samples.theta(:,:,j);
  for c=1:Nchains
    [p, x] = ksdensity(samples.theta(c,:,j));
    plot(x, p, colors(c));
  end
  title(sprintf('theta %d', j));
end
subplot(3,3,9); hold on
for c=1:Nchains
  [p, x] = ksdensity(samples.mu_theta(c,:));
  plot(x, p, colors(c));
end
title(sprintf('mu.theta'))

% Posterior summaries - intervals
figure;
hold on
for j=1:8
  for c=1:Nchains
    q = quantile(samples.theta(c,:,j), [0.1 0.9]);
    h = line([q(1) q(2)], [j j]+c*0.1); set(h, 'color', colors(c));
    q = quantile(samples.theta(c,:,j), [0.5]);
    h=plot(q,j+c*0.1,'*'); set(h, 'color', colors(c));
  end
  legendstr{j} = sprintf('theta %d', j);
end
j = 9;
for c=1:Nchains
  q = quantile(samples.mu_theta(c,:), [0.1 0.9]);
  h = line([q(1) q(2)], [j j]+c*0.1); set(h, 'color', colors(c));
  q = quantile(samples.mu_theta(c,:), [0.5]);
  h=plot(q,j+c*0.1,'*'); set(h, 'color', colors(c));
end
legendstr{j} = sprintf('mu.theta');
set(gca,'yticklabel',[])
xlim = get(gca, 'xlim');
for i=1:j
  text(xlim(1), i, legendstr{i});
end
title('80pc posterior intervals (*=median)')
