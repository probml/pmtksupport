%% Probabilistic Machine Learning (CPSC 530A)
%% Instructor: Nando de Freitas
%% Course Project
%% James Cook, Iryna Skrypnyk, Roland Wenzel

% **********************************
% *** RBPF_Features_to_Obstacles ***
% **********************************

function obs_map = rbpf_features2obstacles(feat_map);

% converts a map of unique features into a map of obstacles

[y,x] = size(feat_map);

for i = 1:y
    for j = 1:x
        if feat_map(i,j) == 0   % this is free space
            obs_map(i,j) = 1;
        else
            obs_map(i,j) = 0;   % this is an obstacle
        end
    end
end
