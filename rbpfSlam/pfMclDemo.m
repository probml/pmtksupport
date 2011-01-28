%% Probabilistic Machine Learning (CPSC 530A)
%% Instructor: Nando de Freitas
%% Course Project
%% James Cook, Iryna Skrypnyk, Roland Wenzel

% ******************************
% *** Particle Filtering     ***
% ******************************

clear all;
close all;
rand('seed', 23);                   % needed for calls to unidrnd
randn('seed', 42);                  % needed for calls to randn
time = cputime;

disp('please wait...')

% SET UP GRAPHICS:
% ===============

set(0,'Units','pixels');
scnsize = get(0,'ScreenSize');
%    [80, 80, ...
figure('Position', ...
  [400, 400, ...
     scnsize(3)-320, ...
     scnsize(4)-320]);              % [left bottom width height]
subplot(1, 1, 1);                   % enable automatic scaling

% LOAD DATA AND SET PARAMETERS:
% ============================

generate_data = 1;                  % flag if data shall be generated again (must be done when changing vision or mode!)
vision = 'f';                       % use vision ('t') or not ('f'); vision requires 2d gaussians
mode = '1d';                        % use 1d ('1d') or 2d ('2d') gaussian to calculate the weights

if (vision == 't')
    mode = '2d';
    disp('using vision.');
else
    disp('using laser scan.');
end
disp(['using ', mode, ' gaussian.']);

N = 300;                            % number of samples (should be higher for vision with a small angle)
disp(['using ', num2str(N), ' particles.']);
z_t = zeros(N,2);                   % samples: z_t(i) = [y(i),x(i)]
w_t = zeros(N,1);                   % weights: w_t(i)

pos_sigma  = 0.3*eye(2);            % Gaussian noise in odometry (used for position proposal)
scan_sigma = 0.2;                   % Gaussian noise in scan readings (will be used twice in 2d case)
angle_sigma = 0.3;                  % Gaussian noise for rotation angle of the robot (only for vision, not used yet)
distance_sigma = 4;                 % used for determining probability of a ray scan

obstacle = 0;                       % obstacle mark in the world of obstacles
ray_length = 10;                    % length of a single ray in each direction
disp(['ray length = ', num2str(ray_length)]);
vision_angle = 90;                  % how much can the robot see?
if (vision == 't')
    disp(['vision angle = ', num2str(vision_angle)]);
end

world = load('pf_world.data');
[WY,WX] = size(world);              % WY(row),WX(col) = dimensions of the world

% generate data if needed (move-, scan-, and positions.data)
if (generate_data == 1)
	disp('generating data...');
	pf_make_data(world, obstacle, ray_length, scan_sigma, mode, vision, vision_angle);
	disp('done');
else
	disp('using old data...');
end

pos = load('pf_positions.data');

delta_pos = load('pf_move.data');   % collection of position changes, one for each move
[M,c] = size(delta_pos);            % M = number of robot moves

scans = load('pf_scan.data');       % collection of ray scans, for each move
if (mode == '2d')
	% segment feature scans for each move
	view = cell(M,1);
    i = 1;  % current entry
    for m = 1:M  % current block
        count = scans(i, 1);
        view{m} = scans(i+1:i+count, :);
        i = i + 1 + count;
    end
    disp(['input consists out of ', num2str(i-M-1), ' scan readings.']);
else
    % no preprocessing needed in 1d mode
    R = 8;                          % number of rays (scan directions)
    view = scans;
end

z_ts = cell(M+1,1);                 % record of position distributions, for each move

robot_color = [0.4, 0.8, 0.3];      % the color of the robot
robot_color_exp = 'red';            % the color of the robot in case of slippage

% distance increments along each ray direction
delta_N  = [-1,  0];
delta_NE = [-1,  1];
delta_E  = [ 0,  1];
delta_SE = [ 1,  1];
delta_S  = [ 1,  0];
delta_SW = [ 1, -1];
delta_W  = [ 0, -1];
delta_NW = [-1, -1];

delta_ray = [delta_N; delta_NE; delta_E; delta_SE; delta_S; delta_SW; delta_W; delta_NW];

% INITIALIZE:
% ==========

for i=1:N                                       % assume initial uniform distribution
    z_t(i,:) = unidrnd([WY,WX]);
    while world(z_t(i,1),z_t(i,2)) == obstacle  % make sure initial proposal is valid
        z_t(i,:) = unidrnd([WY,WX]);
    end
end

z_ts{1} = z_t;

% LET'S GO:
% ========

for m=1:M                                       % for each move
    for i=1:N                                   % for each sample
        z_t(i,:) = ...                          % make a new position proposal
            pos_proposal(z_t(i,:), delta_pos(m,:), world, pos_sigma);
        
    % make a scan from a proposed position
    if (vision == 't')
        alpha = unidrnd(360);                   % "sample" the angle without prior (this could and should be done better...)
       	proposed_view = pf_scan_vision(z_t(i,:), world, ray_length, alpha, vision_angle, scan_sigma);
        disp(['                                                    ', ...
              'checking (', num2str(z_t(i, 2)), '|', num2str(z_t(i, 1)), ') with heading ', num2str(alpha), ' degrees...']);
    else
    	proposed_view = pf_scan_laser(z_t(i,:), world, ray_length, obstacle, scan_sigma, mode);
    end

    if (mode == '1d')
            w_t(i) = ...                        % assign a weight factor to the proposed position
                pf_weight_1d(view(m,:)', proposed_view, distance_sigma);
        else
            w_t(i) = ...                        % assign a weight factor to the proposed position
                pf_weight_2d(view{m}, proposed_view, distance_sigma);
            disp(['weighting move ', num2str(m), ' (', num2str(M), ') / sample ', num2str(i), ...
                  ' (', num2str(N), '): ', num2str(w_t(i))]);
        end
    end

	w_t = w_t./sum(w_t);                        % normalize the weights
	resampled_index = deterministicR(1:N,w_t);
	z_t(:,:) = z_t(resampled_index,:);          % resample  
	z_ts{m+1} = z_t;                            % store current samples
end

disp(['calculations done after ', num2str(cputime-time), ' seconds.']);

% LET'S SEE GRAPHICS:
% ==================

for m=1:(M+1)
    y = pos(m,1);
    x = pos(m,2);

    z_t = z_ts{m};
    %imshow(world);
    imshow(world, 'initialmagnification', 'fit')
    title('Localization with a known map');
    disp_world = repmat(world, [1, 1, 3]);      % get the RGB-values (enables us to use color...)
    for i=1:N                                   % purple = reduce green :-)
        old_color = reshape(disp_world(z_t(i,1), z_t(i,2), :), [1, 3]);
        new_color = max([0, 0, 0], old_color - [0.05, 0.1, 0]);
        disp_world(z_t(i,1), z_t(i,2), :) = new_color;
        rectangle('Position', [z_t(i,2)-0.5, z_t(i,1)-0.5, 1, 1], ...
                  'EdgeColor', new_color, 'FaceColor', new_color);
    end
    
	if (mode == '1d')
        if  (m > 1)
            curr_view = scans(m-1,:);           % display the sensor readings that led to this posterior
            for r = 1:R
                r_x = x + delta_ray(r, 2) * curr_view(r);
                r_y = y + delta_ray(r, 1) * curr_view(r);
                if ((r_x >= 1) & (r_x <= WX) & (r_y >=1) & (r_y <= WY)) % range check; else disp_world will fail
                    rectangle('Position', [r_x-0.5, r_y-0.5, 1, 1], ...
                              'EdgeColor', 'magenta', 'FaceColor', reshape(disp_world(r_y, r_x, :), [1, 3]));
                end
            end
        end
    else
        if  (m > 1)
            curr_view = view{m-1};              % display the sensor readings that led to this posterior
            [V, w] = size(curr_view);
            for v = 1:V
                v_x = x + curr_view(v, 2);
                v_y = y + curr_view(v, 1);
                if ((v_x >= 1) & (v_x <= WX) & (v_y >=1) & (v_y <= WY)) % range check; else disp_world will fail
                    rectangle('Position', [v_x-0.5, v_y-0.5, 1, 1], ...
                              'EdgeColor', 'magenta', 'FaceColor', reshape(disp_world(v_y, v_x, :), [1, 3]));
                end
            end
        end
    end
    
    if (m >= 2)
        % draw a line representing the real movement
        dy = y-old_y;
        dx = x-old_x;
        d = sqrt(dx^2 + dy^2);
        line([old_x, x-0.5/d*dx], [old_y, y-0.5/d*dy], 'Color', robot_color);
    end

    % draw the actual robot position
    rectangle('Position', [x-0.5,      y-0.5,      1,    1],    'Curvature', [1 1], 'EdgeColor', robot_color);
    rectangle('Position', [x-0.5,      y-0.5+0.67, 0.5,  0.5],  'Curvature', [1 1], 'EdgeColor', robot_color);
    rectangle('Position', [x-0.5+0.59, y-0.5+0.84, 0.32, 0.32], 'Curvature', [1 1], 'EdgeColor', robot_color);
    line([x-0.5+0.75, x-0.5+1.00], [y-0.5+0.55, y-0.5+0.60], 'Color', robot_color)
    line([x-0.5+0.75, x-0.5+0.75], [y-0.5+0.30, y-0.5+0.30], 'Color', robot_color)
    
    if (m >= 2)
        % the expected positions
        xx = old_x + delta_pos(m-1, 2);
        yy = old_y + delta_pos(m-1, 1);
	
        if ((yy ~= y) | (xx ~= x)) % slippage!
        	% draw a line representing the odometry reading
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
    end

    Mov(m) = getframe;
    
    old_x = x;
    old_y = y;
    
    drawnow
    pause
end

disp('play movie')
fps=2;
movie(Mov,1,fps);
%movie(Mov, 2, 2);
