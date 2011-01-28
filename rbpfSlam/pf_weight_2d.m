%% Probabilistic Machine Learning (CPSC 530A)
%% Instructor: Nando de Freitas
%% Course Project
%% James Cook, Iryna Skrypnyk, Roland Wenzel

% ***********************************
% *** PF_Weight (for 2d gaussian) ***
% ***********************************

function w = pf_weight_2d(curr_view, proposed_view, dis_sigma)

w = 1;
[D,r] = size(curr_view);
[E,r] = size(proposed_view);

% compute joint probability (assume independency)
for i=1:D   % for each obstacle in the current view (don't check those expected ones that you don't see)
    closest = 1;
    min_distance = 1e99;
    
    for j=1:E   % find the closest one in the proposition
        delta = curr_view(i,:) - proposed_view(j,:);
        distance = sqrt(delta(1)^2 + delta(2)^2);
        
        if (distance < min_distance)
            closest = j;
            min_distance = distance;
        end
    end
    
    % and include the relevant weight
    w = w * gaussian(min_distance, dis_sigma);  % "independent"!
end


% ******************************
% *** Gaussian               ***
% ******************************

function p = gaussian(x_minus_mu, sigma)

% Gaussian distribution (up to a normalization constant)
p = exp(-1/(2*sigma^2)*(x_minus_mu)^2) + 1e-99;
