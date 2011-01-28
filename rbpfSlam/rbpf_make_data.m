%% Probabilistic Machine Learning (CPSC 530A)
%% Instructor: Nando de Freitas
%% Course Project
%% James Cook, Iryna Skrypnyk, Roland Wenzel

% ******************************
% *** RBPF_Make_Data         ***
% ******************************

function rbpf_make_data(obs_world, free_space, ray_length, scan_sigma)

[feat_world, n_features] = rbpf_obstacles2features(obs_world);

positions = ...                     % robot path for localization and mapping
    [24,15;                         % initial location (which the robot knows)
     22,15; 21,14; 19,14; 19,12; 21,12; 22,13; 23,13; 23,11; 24, 9; 23, 7;
     21, 7; 19, 7; 17, 7; 16, 9; 17,11; 17,13; 16,14; 14,14; 13,14; 13,12;
     13,10; 14, 8; 13, 7; 12, 8; 12,10; 10,10;  9, 9;  9, 7;  9, 5;  8, 6;
      8, 8;  8,10;  9,11;  9,13;  8,14;  6,14;  6,12;  6,10;  8,12; 10,12;
     12,11; 14,11; 15, 9; 16, 7; 16, 5; 18, 5; 20, 6; 22, 5; 24, 5; 24, 7;
     23, 9; 24,10; 24,12; 24,13; 24,15];

[r,c] = size(positions);

deltas = [];
scans  = [];
for i=2:r                           % generate move deltas and scans along the path
	delta = positions(i,:) - positions(i-1,:);
	deltas = [deltas; delta];

	scan = rbpf_scan(positions(i,:), feat_world, ray_length, free_space, scan_sigma);
	scans = [scans; scan];
end
  
% noise in odometry
deltas( 3,2) = deltas( 3,2)+1;
deltas( 4,1) = deltas( 4,1)-3;
deltas( 4,2) = deltas( 4,2)-3;
deltas(12,1) = deltas(12,1)-1;
deltas(12,2) = deltas(12,2)+1;
deltas(21,1) = deltas(21,1)-1;
deltas(24,2) = deltas(24,2)+1;
deltas(33,1) = deltas(33,1)+1;
deltas(41,2) = deltas(41,2)-1;
deltas(45,1) = deltas(45,1)-1;
deltas(48,1) = deltas(48,1)+1;

save('rbpf_positions.data', 'positions', '-ASCII');
save('rbpf_move.data',      'deltas',    '-ASCII');
save('rbpf_scan.data',      'scans',     '-ASCII');
