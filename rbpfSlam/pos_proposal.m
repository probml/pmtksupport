%% Probabilistic Machine Learning (CPSC 530A)
%% Instructor: Nando de Freitas
%% Course Project
%% James Cook, Iryna Skrypnyk, Roland Wenzel

% ******************************
% *** Pos_Proposal           ***
% ******************************

function new_pos = pos_proposal(old_pos, delta_pos, world, pos_sigma)

desired_pos = round(old_pos + delta_pos + randn(size(old_pos)) * pos_sigma);

world_size = size(world);

% The following code ensures the robot does not move through obstacles.
% The classic differential line drawing algorithm with modifications.

dx = desired_pos(2) - old_pos(2);
dy = desired_pos(1) - old_pos(1);

ax = abs(dx)*2;
ay = abs(dy)*2;

sx = sign(dx);
sy = sign(dy);

% Coordinates drawn out by the differential algorithm.
% These coordinates can pass through walls. 
tracex = old_pos(2);
tracey = old_pos(1);

% The actual coordinates. We don't increase these if the robot will run into
% a wall.
actualx = old_pos(2);
actualy = old_pos(1);

% check if we are on the edge of the world
if (actualx + sx < 1) | (actualx + sx > world_size(2)) ...
    | (actualy + sy < 1) | (actualy + sy > world_size(1))
	% if so, get out
	new_pos = old_pos;
	%break; 
	return;
end

if (ax > ay)   % x dominant. moving mostly in x direction
	d = ay-floor(ax/2);
	
	while (tracex ~= desired_pos(2))
		if (d >= 0)
			tracey = tracey + sy;
			d = d - ax;
			otherstep = 1; 
		else
			otherstep = 0;
		end
		tracex = tracex + sx;
		d = d + ay;
		
		% Now we check if the incremental step will put the robot inside an obstacle
		if world(actualy + otherstep*sy, actualx + sx) == 0
			if otherstep == 1
				% We are moving diagonally and there is an obstacle on the diagonal.  
				if (world(actualy, actualx + sx) ~= 0)
					% There is an obstacle in the diagonal but not on the horizontal 
					% Move horizontally instead
					actualx = actualx + sx;
				elseif (world(actualy + sy, actualx) ~= 0)
					% There is an obstacle in the diagonal and horizontal, but not vertical
					% Move vertically instead
					actualy = actualy + sy;
				end
		    end 
		else
			% else we are moving only horizontally and
			%      we hit an obstacle. we don't increase movement
	
            % No obstacle in the desired step 
			actualy = actualy + otherstep*sy;
			actualx = actualx + sx;
		end
    end
    
else % y dominant. moving mostly in y direction
	d = ax - floor(ay/2); 
	while (tracey ~= desired_pos(1))
		if (d >= 0)
			tracex = tracex + sx;
			d = d - ay;
			otherstep = 1;
		else
			otherstep = 0;
		end
		tracey = tracey + sy;
		d = d + ax;
		
		% Now we check if the incremental step will put the robot inside an obstacle
		if world(actualy + sy, actualx + otherstep*sx) == 0
			if otherstep == 1
			    % We are moving in a diagonal direction and there is an obstacle on the diagonal.
				if (world(actualy + sy, actualx) ~= 0)
					% There is an obstacle in the diagonal but not on the vertical
					% Move vertically instead
					actualy = actualy + sy;
				elseif (world(actualy, actualx + sx) ~= 0) 
					% There is an obstacle in the diagonal and vertical, but not horizontal 
					% Move horizontally instead
					actualx = actualx + sx;
				end
			end
		else
		% else we are moving only in the vertical direction and
		%      we hit an obstacle. don't increase movement
	
        % No obstacle in the desired step
			actualy = actualy + sy;
			actualx = actualx + otherstep*sx;
		end
	end
end

new_pos = [actualy actualx];
