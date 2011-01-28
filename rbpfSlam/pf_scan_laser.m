%% Probabilistic Machine Learning (CPSC 530A)
%% Instructor: Nando de Freitas
%% Course Project
%% James Cook, Iryna Skrypnyk, Roland Wenzel

% ******************************
% *** PF_Scan_Laser          ***
% ******************************

function view = pf_scan_laser(curr_pos, world, ray_length, obstacle, scan_sigma, mode)

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

[D,c] = size(delta_ray);                % number of individual rays (scan directions)

distances = ones(1,D)*1e99;             % distances to obstacles in each ray direction
                                        % default distance (no obstacle detected) = 1e99

obstacles = zeros(1,D);                 % obstacle detection flags

ray_pos = zeros(D,2);                   % current ray extent in each direction
ray_pos(:,1) = curr_pos(1);
ray_pos(:,2) = curr_pos(2);

for r=1:ray_length                      % at the current ray extent     
    for d=1:D                           % for each direction
        if obstacles(d) == 0            % no obstacle detected yet in this direction
            ray_pos(d,:) = ray_pos(d,:) + delta_ray(d,:);
                        
            y = ray_pos(d,1);
            x = ray_pos(d,2);
            
            if world(y,x) == obstacle   % obstacle detected in the world
                obstacles(d) = 1;
                distances(d) = r;
            end
        end
    end
end

if (mode == '1d')
    % use 1d gaussian
    view = distances;
    % add 1d noise
    % those directions where we can't see anything in the laser range will be ignored in the
    % 1 dimensional weight function
    for d = 1:D
        distances(d) = add_scan_noise_1d(distances(d), ray_length, scan_sigma);
    end
else
    % use 2d gaussian
    view = [];
    for d=1:D                           % return a vector of detected obstacles rather than just distances
        % rebuild and add 2d noise
        % ignore those directions where we can't see anything in the laser range
        if (distances(d) ~= 1e99)
            view = [view; add_scan_noise_2d(delta_ray(d, :) * distances(d), ray_length, scan_sigma)];
        end
    end
end


% ******************************
% *** Add_Scan_Noise         ***
% ******************************

function new_distance = add_scan_noise_1d(old_distance, ray_length, scan_sigma)

new_distance = round(old_distance + scan_sigma * randn(1,1));
    
if new_distance < 1                     % stay within boundaries
    new_distance = 1;        
end
if new_distance > ray_length            % stay within boundaries
    new_distance = ray_length;
end


function new_position = add_scan_noise_2d(old_position, ray_length, scan_sigma)

new_x = round(old_position(2) + scan_sigma * randn(1,1));
new_y = round(old_position(1) + scan_sigma * randn(1,1));
    
% if new_x < 1                            % stay within boundaries
%     new_x = 1;        
% end    
% if new_x > ray_length                   % stay within boundaries
%     new_x = ray_length;
% end
% 
% if new_y < 1                            % stay within boundaries
%     new_y = 1;        
% end    
% if new_y > ray_length                   % stay within boundaries
%     new_y = ray_length;
% end

new_position = [new_y, new_x];
