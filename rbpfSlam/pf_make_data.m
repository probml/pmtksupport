%% Probabilistic Machine Learning (CPSC 530A)
%% Instructor: Nando de Freitas
%% Course Project
%% James Cook, Iryna Skrypnyk, Roland Wenzel

% ******************************
% *** PF_Make_Data           ***
% ******************************

function pf_make_data(world, obstacle, ray_length, scan_sigma, mode, vision, vision_angle)

positions = ...                     % robot path for localization
    [16, 7;                         % initial location (of which the robot is unaware)
     14, 7; 12, 7; 12, 9; 12,11; 12,14; 10,14;  8,14; 5,13; 4,16; 6,17;
      8,17; 10,17; 10,20; 11,22;  9,24;  8,26;  8,29];
 
[r,c] = size(positions);

deltas = [];
scans  = [];
for i=2:r                           % generate deltas and scans along the path
    delta = positions(i,:) - positions(i-1,:);
	deltas = [deltas; delta];

    if (vision == 't')
        if (delta(1) < 0)           % robot goes to the north
            alpha = acos(delta(2)/sqrt(delta(1)^2+delta(2)^2)) /pi*180;
        else
            alpha = 360 - acos(delta(2)/sqrt(delta(1)^2+delta(2)^2)) /pi*180;
        end
    	scan = pf_scan_vision(positions(i,:), world, ray_length, alpha, vision_angle, scan_sigma);
    else
    	scan = pf_scan_laser(positions(i,:), world, ray_length, obstacle, scan_sigma, mode);
    end
    if (mode == '2d')               % arbitrary number of features in 2d mode
    	scans = [scans; size(scan, 1), -1; scan];
    else                            % always use 8 rays in 1d mode
    	scans = [scans; scan];
    end
end

% noise in odometry
deltas( 2,2) = deltas( 2,2)+1;
deltas( 6,1) = deltas( 6,1)+1;
deltas(11,1) = deltas(11,1)+1;
deltas(11,2) = deltas(11,2)-1;
deltas(16,1) = deltas(16,1)+1;

save('pf_positions.data', 'positions', '-ASCII');
save('pf_move.data',      'deltas',    '-ASCII');
save('pf_scan.data',      'scans',     '-ASCII');
