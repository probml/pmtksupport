%% Probabilistic Machine Learning (CPSC 530A)
%% Instructor: Nando de Freitas
%% Course Project
%% James Cook, Iryna Skrypnyk, Roland Wenzel

% ******************************
% *** PF_scan_vision         ***
% ******************************

function curr_view = pf_scan_vision(position, world, range, heading, degrees, scan_sigma)

draw = 0;
check_occluded = 1;

xpos = position(2);
ypos = position(1);

[ysize, xsize] = size(world);

M = [];                                     % no known object to far...
world_in_range = ones (ysize, xsize) / 2;   % grey
beams = ones (ysize, xsize) / 2;            % grey (again)
checked = ones (ysize, xsize) / 2;          % grey (and again)
view = ones (ysize, xsize) / 2;             % grey (still)

for j = 1:ysize
    for i = 1:xsize
        [alpha, d] = cart2pol(i-xpos, j-ypos);                              % relative to the current position
        rel_0based_degree = mod(alpha/pi*180 + heading + degrees/2, 360);     % 0..(degrees) are the perceived angles
        if (rel_0based_degree < degrees) & (d < range)                      % it is in range and in a good angle
            % remember this coordinate in the matrix (obstacle or non-obstacle, occluded or not occluded)
            world_in_range(j, i) = world(j, i);
            M = [M; j-ypos, i-xpos, alpha, d, world(j, i)];
        end
    end
end
% mark the robots position
world_in_range(ypos, xpos) = 0.8;
checked(ypos, xpos) = 0.8;
view(ypos, xpos) = 0.8;

if (draw == 1)
	set(0,'Units','pixels') 
	scnsize = get(0,'ScreenSize');
	figure('Position', [80, 80, scnsize(3) - 160, scnsize(4) - 160]); % [left bottom width height]

    subplot(2, 3, 1)    % world
	imshow (world)
    title('world')
    subplot(2, 3, 2)    % world_in_range
	imshow (world_in_range)
    title('world in range')
    subplot(2, 3, 3)    % info
    imshow(ones(ysize, xsize) * 0.7);
    text(3, 3, ['X = ', num2str(xpos)]);
    text(3, 5, ['Y = ', num2str(ypos)]);
    title('info')
end

M = sortrows(M, 4);     % sort by distance (closer objects occlude those that are farther away)
                        % do sorting only once, since the objects don't move...
[r, c] = size(M);
for j = 1:r             % iterate through all points
    y     = M(j, 1);    % the obstacle in the "beam" in relative cartesian coordinates
    x     = M(j, 2);
    alpha = M(j, 3);    % angle to the occluding/reference point
    d     = M(j, 4);    % distance to the point
    color = M(j, 5);    % "color" (object/free, occluded/not occluded) of that point
    if d ~= 0           % this point is not the where the robot stands
        if (color == 0) | ((check_occluded == 1) & (color == 0.3))   % this is an obstacle and occludes other points; mask all fields behind this one

            if (draw == 1)
                subplot(2, 3, 4)    % beams
				imshow (beams)
                title('beams')
                subplot(2, 3, 5)    % checked
				imshow (checked)
                title('checked points')
                subplot(2, 3, 6)    % view
				imshow (view)
                title('perception')

                subplot(2, 3, 4)    % beams
                for i = 1:j         % those points are closer than the occluding point; just draw them
                    b_y  = M(i, 1); % cartesian coordinates of the eventually occluded point
                    b_x  = M(i, 2);
                    H = rectangle('Position', [xpos+b_x-0.5, ypos+b_y-0.5, 1, 1]);
                    set(H, 'EdgeColor', 'white', 'FaceColor', ones(1, 3) * M(i, 5));
                end
            end

            % find the relevant corners of the occluding "box"
            % this is important to do it "right"! however, it can be approximated much faster, but will cause
            % misperceptions then...
            % r_1 should be the point with the smaller angle (clockwise (!) spoken, because the y's increase to the south)
            if (y > 0.5)
                if (x > 0.5)        % 1
                    %disp('case 1');
                    r_1 = [x+0.5, y-0.5];
                    r_2 = [x-0.5, y+0.5];
                elseif (x >= -0.5)  % 2
                    %disp('case 2');
                    r_1 = [x+0.5, y-0.5];
                    r_2 = [x-0.5, y-0.5];
                else                % 3
                    %disp('case 3');
                    r_1 = [x+0.5, y+0.5];
                    r_2 = [x-0.5, y-0.5];
                end
            elseif (y >= -0.5)
                if (x > 0.5)        % 4
                    %disp('case 4');
                    r_1 = [x-0.5, y-0.5];
                    r_2 = [x-0.5, y+0.5];
                elseif (x >= -0.5)  % 5
                    % should not happen... that's me!
                    % normally I would need 4 point in that case and everything would be occluded
                    %disp('case 5');
                    r_1 = [x-0.5, y+0.5];
                    r_2 = [x+0.5, y-0.5];
                else                % 6
                    %disp('case 6');
                    r_1 = [x+0.5, y+0.5];
                    r_2 = [x+0.5, y-0.5];
                end
            else
                if (x > 0.5)        % 7
                    %disp('case 7');
                    r_1 = [x-0.5, y-0.5];
                    r_2 = [x+0.5, y+0.5];
                elseif (x >= -0.5)  % 8
                    %disp('case 8');
                    r_1 = [x-0.5, y+0.5];
                    r_2 = [x+0.5, y+0.5];
                else                % 9
                    %disp('case 9');
                    r_1 = [x-0.5, y+0.5];
                    r_2 = [x+0.5, y-0.5];
                end
            end

            % now calculate the minimum angle covering the whole reference point (seen from (0|0))
            [alpha_1, a1_d] = cart2pol(r_1(1), r_1(2));
            [alpha_2, a2_d] = cart2pol(r_2(1), r_2(2));

            if (draw == 1)
                subplot(2, 3, 4)    % beams
                % the robot
   				H = rectangle('Position', [xpos-0.5, ypos-0.5, 1, 1]);
                set(H, 'EdgeColor', 'red', 'Curvature', [1 1]);
                % the occluding point
				H = rectangle('Position', [xpos+x-0.5, ypos+y-0.5, 1, 1]);
                set(H, 'EdgeColor', 'green');
                % the left beam
                H = line([xpos, xpos+r_1(1)], [ypos, ypos+r_1(2)]);
                set(H, 'Color', 'green');
                % the right beam
                H = line([xpos, xpos+r_2(1)], [ypos, ypos+r_2(2)]);
                set(H, 'Color', 'blue');
            end
            
            checked(ypos+y, xpos+x) = color;
            if (draw == 1)
                subplot(2, 3, 5)    % checked
                H = rectangle('Position', [xpos+x-0.5, ypos+y-0.5, 1, 1]);
                set(H, 'EdgeColor', 'red', 'FaceColor', ones(1, 3) * color);
            end
            if (color == 0) % black
                view(ypos+y, xpos+x) = color;
                if (draw == 1)
                    subplot(2, 3, 6)    % view
                    H = rectangle('Position', [xpos+x-0.5, ypos+y-0.5, 1, 1]);
                    set(H, 'EdgeColor', 'red', 'FaceColor', ones(1, 3) * color);
                end
            end

            for i = j+1:r
                b_y  = M(i, 1);     % cartesian coordinates of the eventually occluded point
                b_x  = M(i, 2);
                beta = M(i, 3);     % angle to the (occluded/not occluded) point

                if (beta > alpha_1) & (beta < alpha_2)  % occluded; Warning: this line might be buggy!
                    if M(i, 5) == 0
                        M(i, 5) = 0.3;
                    elseif M(i, 5) == 1
                        M(i, 5) = 0.7;
                    % else: already occluded
                    end
                end

                if (draw == 1)
                    subplot(2, 3, 4)    % beams
                    % the eventually occluded point
                    H = rectangle('Position', [xpos+b_x-0.5, ypos+b_y-0.5, 1, 1]);
                    set(H, 'EdgeColor', 'yellow', 'FaceColor', ones(1, 3) * M(i, 5));
                end
                
            end % for all farther points than the occluding one

            if (draw == 1)
                drawnow;
            end
            
        end % if there is an obstacle at this point (that can occlude other points)
    end % if the point is not the robot
end % iterate through all points

if (draw == 1)
	subplot(2, 3, 5)    % checked
	imshow (checked)
	title('checked points')
	subplot(2, 3, 6)    % view
	imshow (view)
	title('perception')
end

% generate the list of obstacles
curr_view = [];

for j = 1:ysize
    for i = 1:xsize
        if (view(j, i) == 0)
            curr_view = [curr_view; add_scan_noise_2d([j-ypos, i-xpos], range, scan_sigma)];
        end
    end
end


% ******************************
% *** Add_Scan_Noise         ***
% ******************************

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
