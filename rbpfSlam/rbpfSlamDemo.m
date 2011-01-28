%% Probabilistic Machine Learning (CPSC 530A)
%% Instructor: Nando de Freitas
%% Course Project
%% James Cook, Iryna Skrypnyk, Roland Wenzel

% ********************************************
% *** Rao-Blackwellised Particle Filtering ***
% ********************************************

clear all;
close all;
rand('seed', 23);                   % needed for calls to unidrnd
randn('seed', 42);                  % needed for calls to randn

disp('please wait...');

% SET UP GRAPHICS:
% ===============

set(0,'Units','pixels');
scnsize = get(0,'ScreenSize');
figure('Position', ...
    [80, 80, ...
     scnsize(3)/2-200, ...
     scnsize(4)-320]);              % [left bottom width height]
subplot(1, 1, 1);                   % enable automatic scaling for figure 1
figure('Position', ...
    [scnsize(3)/2+80, 80, ...
     scnsize(3)/2-200, ...
     scnsize(4)-320]);              % [left bottom width height]
subplot(1, 1, 1);                   % enable automatic scaling for figure 2


% LOAD DATA AND SET PARAMETERS:
% ============================

generate_data = 1;                  % flag if data shall be generated again

N = 50;                            % number of samples
z_t = zeros(N,2);                   % location samples: z_t(i) = [y(i),x(i)]
x_t = cell(N,1);                    % map samples:      x_t{i} = [detected(j),y(j),x(j)]^F; F = number of features
w_t = zeros(N,1);                   % weights:          w_t(i)

pos_sigma  = 0.5*eye(2);            % Gaussian noise in odometry (y, x)
scan_sigma = 0.2*eye(2);            % Gaussian noise in scan readings (y, x)

free_space = 0;                     % free space mark in the world of *features*
ray_length = 10;                    % length of a single ray in each direction

obs_world = load('rbpf_world.data');                    % load the world
[feat_world,F] = rbpf_obstacles2features(obs_world);    % F = number of features

% generate data if needed (move-, scan-, and positions.data)
if (generate_data == 1)
	disp('generating data...');
	rbpf_make_data(obs_world, free_space, ray_length, scan_sigma);
	disp('done');
end

delta_pos = load('rbpf_move.data'); % collection of position changes, one for each move
[M,c] = size(delta_pos);            % M = number of robot moves

R = 8;                              % number of rays (scan directions)
scans = load('rbpf_scan.data');     % collections of features and their relative distances
% segment feature scans for each move
view = cell(M,1);
for i=1:M
    view{i} = scans(R*i-R+1:R*i,:); % i-th block of R data each
end

% Initialize parameters for movies
pos          = load('rbpf_positions.data');  
init_pos     = pos(1,:);
old_x        = init_pos(2);
old_y        = init_pos(1);
[WY,WX]      = size(obs_world);

robot_color_real = [0.8, 0.3, 1];   % the color of the robot in figure 1
robot_color_view = [0.3, 0.6, 0.2]; % the color of the robot in figure 2
robot_color_exp  = 'red';           % the color of the robot in figure 2 in case of slippage

% INITIALIZE:
% ==========

init_view = rbpf_scan(init_pos, feat_world, ray_length, free_space, scan_sigma);

B = 0.3*eye(F);
D = 0.2*eye(F);                     % = scan_sigma

for i=1:N                           % starting at a known position with initial view
	z_t(i,:)     = init_pos;
	x_t{i}       = zeros(F,3);      % first element = 0: feature not detected yet
	x_t{i}       = rbpf_add_features(x_t{i}, init_view);
	sigma(:,:,i) = 1*eye(F);
end  

% LET'S GO:
% ========

% remark: the world before the first move will not be drawn.

for m=1:M                                       % for each move
  for i=1:N                                   % for each sample
    old_pos = z_t(i,:);
    
    z_t(i,:) = ...                          % make a new position proposal
	pos_proposal(z_t(i,:), delta_pos(m,:), obs_world, pos_sigma);
    
    x_t{i} = ...                            % update the relative position of the existing features
	rbpf_update_features(x_t{i}, z_t(i,:) - old_pos);
    
    % KALMAN PREDICTION: (from Nando's demo_rbpf_gauss.m)
    % =================
    
    % Our C selects features from the map that the robot currently sees.
    % The observation model is y = C*x_t + N(0,D) (y is in the map format).
    % C selects features from x_t to get y. If feature i is seen in view{m}
    % then we put a 1 on the diagonal at C(i,i). The rest of C is zero.
    
    C = zeros(F);                           % selector
    for r=1:R                               % for each ray (scan direction)
      f = view{m}(r,1);                   % feature id
      if (f ~= 0) & (x_t{i}(f,1) == 1)    % feature has been detected before
	C(f,f) = 1;                     %   then: select it
      end;
    end
    
    sigma(:,:,i)  = sigma(:,:,i)+C*B*B';
    St(:,:,i)     = C*sigma(:,:,i)*C'+D*D';
    y_pred(:,:,i) = C*x_t{i}(:,2:3);        % observed features from the map
    
  end
  
  y_view = zeros(F,2);                        % observed features from the current view
  for r=1:R                                   % for each ray (scan direction)
    f = view{m}(r,1);                       % feature id
    if (f ~= 0)
      y_view(f,1) = view{m}(r,2);
      y_view(f,2) = view{m}(r,3); 
    end
    % all unknown features are still (0, 0) (which is invalid), so those entries
    % will be ignored in the weight-function
  end  
  
  % CALCULATE THE WEIGHTS:
  % =====================
  
  for i = 1:N  
    w_t(i) = rbpf_weight(y_view, y_pred(:,:,i), St(:,:,i));
  end
  w_t = w_t./sum(w_t);                        % normalize the weights
  
  % RESAMPLE:
  % ========
  
  resampled_index = deterministicR(1:N,w_t);
  
  z_t(:,:) = z_t(resampled_index,:);        
  x_t      = x_t(resampled_index);               
  sigma    = sigma(:,:,resampled_index);
  St       = St(:,:,resampled_index);
  y_pred   = y_pred(:,:,resampled_index); 
  
  % UPDATE THE MAP: (with Kalman filter)
  % ==============
  
  for i=1:N
    K             = sigma(:,:,i)*C'*pinv(St(:,:,i));
    x_t{i}(:,2:3) = x_t{i}(:,2:3)+K*(y_view-y_pred(:,:,i));
    sigma(:,:,i)  = sigma(:,:,i)-K*C*sigma(:,:,i);
  end
  
  % ADD NEW FEATURES TO THE MAP:
  % ===========================
  
  for i=1:N                                   % for each surviving sample
    x_t{i} = rbpf_add_features(x_t{i}, view{m});    % known ones will be ignored
  end
  
  % LET'S SEE MOVIES:
  % ================
  
  % Figure 1 displays the map and the average of the locations.
  % It also displays the actual robot position (in robot_color_real)).
  figure(1);
  
  disp_world = obs_world;
  
  for j=1:N   % plot the particles
    disp_world(z_t(j,1),z_t(j,2)) = disp_world(z_t(j,1),z_t(j,2)) - w_t(j);
  end
  
  disp_world = floor(disp_world.*256);
  %imshow(disp_world,gray(256));
  imshow(disp_world,gray(256),'initialmagnification','fit');
  title('What is out there');
  
  % the new (real) positions
  x = pos(m+1,2);
  y = pos(m+1,1);
  
  % draw the actual robot position
  rectangle('Position', [x-0.5,      y-0.5,      1,    1],    'Curvature', [1 1], 'EdgeColor', robot_color_real);
  rectangle('Position', [x-0.5,      y-0.5+0.67, 0.5,  0.5],  'Curvature', [1 1], 'EdgeColor', robot_color_real);
  rectangle('Position', [x-0.5+0.59, y-0.5+0.84, 0.32, 0.32], 'Curvature', [1 1], 'EdgeColor', robot_color_real);
  line([x-0.5+0.75, x-0.5+1.00], [y-0.5+0.55, y-0.5+0.60], 'Color', robot_color_real);
  line([x-0.5+0.75, x-0.5+0.75], [y-0.5+0.30, y-0.5+0.30], 'Color', robot_color_real);
  Mov(m) = getframe;
  
  % Figure 2 displays the map distribution, the real
  % robot position (in robot_color_view) and
  % the expected position of the robot measured by
  % odometry (in robot_color_exp).
  figure(2);
  
  robot_map_avg = zeros(WY,WX);
  
  for j = 1:N
    current_map = x_t{j};
    for k = 1:F
      if x_t{2}(k) ~= 0                                   % if we have seen this feature
	x = round(current_map(k,3)) + z_t(j,2);         % Get absolute positions
	y = round(current_map(k,2)) + z_t(j,1);
	if (y >= 1) & (y <= WY) & (x >=1) & (x <= WX)   % Check that we are not outside bounds
	  robot_map_avg(y,x) = robot_map_avg(y,x) + w_t(j); 
	end
      end
    end
  end
  
  robot_map_avg = ones(size(robot_map_avg)) - robot_map_avg / max(max(robot_map_avg));
  robot_map_avg = floor(robot_map_avg.* 256);
  imshow(robot_map_avg,gray(256),'initialmagnification','fit');
  title('What the robot sees');
  
  % the real position
  x = pos(m+1,2);
  y = pos(m+1,1);
  
  % draw a line representing the real movement
  dy = y-old_y;
  dx = x-old_x;
  d = sqrt(dx^2 + dy^2);
  line([old_x, x-0.5/d*dx], [old_y, y-0.5/d*dy], 'Color', robot_color_view);
  
  % draw the actual robot position
  rectangle('Position', [x-0.5,      y-0.5,      1,    1],    'Curvature', [1 1], 'EdgeColor', robot_color_view);
  rectangle('Position', [x-0.5,      y-0.5+0.67, 0.5,  0.5],  'Curvature', [1 1], 'EdgeColor', robot_color_view);
  rectangle('Position', [x-0.5+0.59, y-0.5+0.84, 0.32, 0.32], 'Curvature', [1 1], 'EdgeColor', robot_color_view);
  line([x-0.5+0.75, x-0.5+1.00], [y-0.5+0.55, y-0.5+0.60], 'Color', robot_color_view);
  line([x-0.5+0.75, x-0.5+0.75], [y-0.5+0.30, y-0.5+0.30], 'Color', robot_color_view);
  
  % the expected positions
  xx = old_x + delta_pos(m, 2);
  yy = old_y + delta_pos(m, 1);
  
  if ((yy ~= y) | (xx ~= x)) % slippage!     draw a line representing the odometry reading
    dy = yy-old_y;
    dx = xx-old_x;
    d = sqrt(dx^2 + dy^2);
    line([old_x, xx-0.5/d*dx], [old_y, yy-0.5/d*dy], 'Color', robot_color_exp);
    
    % draw the expected robot position
    rectangle('Position', [xx-0.5,      yy-0.5,      1,    1],    'Curvature', [1 1], 'EdgeColor', robot_color_exp);
    rectangle('Position', [xx-0.5,      yy-0.5+0.67, 0.5,  0.5],  'Curvature', [1 1], 'EdgeColor', robot_color_exp);
    rectangle('Position', [xx-0.5+0.59, yy-0.5+0.84, 0.32, 0.32], 'Curvature', [1 1], 'EdgeColor', robot_color_exp);
    line([xx-0.5+0.75, xx-0.5+1.00], [yy-0.5+0.55, yy-0.5+0.60], 'Color', robot_color_exp);
    line([xx-0.5+0.75, xx-0.5+0.75], [yy-0.5+0.30, yy-0.5+0.30], 'Color', robot_color_exp);
  end
  
  % plot the observed obstacles
  for i = 1:R % for all rays
    loc = view{m}(i,:,:);
    if ((loc(2) ~= 0) | (loc(3) ~= 0))
      rectangle('Position', [xx+loc(3)-0.5, yy+loc(2)-0.5, 1, 1], ...
		'EdgeColor', 'magenta', ...
		'FaceColor', ones(1, 3) * robot_map_avg(yy+loc(2), xx+loc(3))/256);
    end
  end
  
  Mov2(m) = getframe;
  
  old_x = x;
  old_y = y;
  
  drawnow
  pause
end
