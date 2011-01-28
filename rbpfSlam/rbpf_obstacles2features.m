%% Probabilistic Machine Learning (CPSC 530A)
%% Instructor: Nando de Freitas
%% Course Project
%% James Cook, Iryna Skrypnyk, Roland Wenzel

% **********************************
% *** RBPF_Obstacles_to_Features ***
% **********************************

function [feat_map, n_features] = rbpf_obstacles2features(obs_map);

% converts a map of obstacles into a map of unique features

[y,x] = size(obs_map);
next_unique_number = 1;

for i = 1:y
    for j = 1:x
        if obs_map(i,j) == 0    % this is an obstacle
            feat_map(i,j)      = next_unique_number;
            next_unique_number = next_unique_number + 1;
        else
            feat_map(i,j) = 0;  % this is free space
        end
    end
end

n_features = next_unique_number-1;
