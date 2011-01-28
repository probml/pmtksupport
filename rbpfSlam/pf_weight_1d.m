%% Probabilistic Machine Learning (CPSC 530A)
%% Instructor: Nando de Freitas
%% Course Project
%% James Cook, Iryna Skrypnyk, Roland Wenzel

% ***********************************
% *** PF_Weight (for 1d gaussian) ***
% ***********************************

function w = pf_weight_1d(curr_view, proposed_view, sigma)

w = 1;
[D,c] = size(curr_view);

% compute joint probability (assume independency)
for d = 1:D             % for each ray (scan direction)
    if (curr_view(d) < 1e90 & proposed_view(d) < 1e90)  % there is noise on the marker
        w = w * gaussian(curr_view(d), proposed_view(d), sigma);
    else                % don't include this direction as we can't compare the readings
        w = w * 0.8;    % penalise this reading
    end
end
    
% ******************************
% *** Gaussian               ***
% ******************************

function p = gaussian(x, mu, sigma)

% Gaussian distribution (up to a normalization constant)
p = exp(-1/(2*sigma^2)*(x - mu)^2) + 1e-99;
