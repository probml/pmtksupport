%% Probabilistic Machine Learning (CPSC 530A)
%% Instructor: Nando de Freitas
%% Course Project
%% James Cook, Iryna Skrypnyk, Roland Wenzel

% ******************************
% *** RBPF_Scan              ***
% ******************************

function view = rbpf_scan(curr_pos, world, ray_length, free_space, scan_sigma)

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

distances = zeros(D,3);                 % distances to features in each ray direction

features = zeros(1,D);                  % feature detection flags

ray_pos = zeros(D,2);                   % current ray extent in each direction
ray_pos(:,1) = curr_pos(1);
ray_pos(:,2) = curr_pos(2);

for r=1:ray_length                      % at the current ray extent     
    for d=1:D                           % for each direction
        if features(d) == 0             % no feature detected yet in this direction
            ray_pos(d,:) = ray_pos(d,:) + delta_ray(d,:);
                        
            y = ray_pos(d,1);
            x = ray_pos(d,2);
            
            if world(y,x) ~= free_space % feature detected in the world
                features(d) = 1;
                distances(d,:) = [world(y,x),y-curr_pos(1),x-curr_pos(2)];
                distances(d,2:3) = add_scan_noise(distances(d,2:3), ray_length, scan_sigma);
            end
        end
    end
end

view = distances;

% ******************************
% *** Add_Scan_Noise         ***
% ******************************

function new_distance = add_scan_noise(old_distance, ray_length, scan_sigma)

new_distance = round(old_distance + randn(size(old_distance)) * scan_sigma);
    
for i=1:2
    if new_distance(i) > ray_length         % stay within boundaries
        new_distance(i) = ray_length;
    elseif new_distance(i) < -ray_length   % stay within boundaries
        new_distance(i) = -ray_length;
    end
end
