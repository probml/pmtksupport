%% Probabilistic Machine Learning (CPSC 530A)
%% Instructor: Nando de Freitas
%% Course Project
%% James Cook, Iryna Skrypnyk, Roland Wenzel

% ******************************
% *** RBPF_Update_Features   ***
% ******************************

function new_map = rbpf_update_features(old_map, delta_pos)

% updates existing features in the current map

new_map = old_map;

[F,c] = size(new_map);

for f=1:F                       % for each possible feature
    if new_map(f,1) == 1        % feature has been detected before
        new_map(f,2:3) = new_map(f,2:3) - delta_pos;
    end
end
