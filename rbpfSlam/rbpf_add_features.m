%% Probabilistic Machine Learning (CPSC 530A)
%% Instructor: Nando de Freitas
%% Course Project
%% James Cook, Iryna Skrypnyk, Roland Wenzel

% ******************************
% *** RBPF_Add_Features      ***
% ******************************

function new_map = rbpf_add_features(old_map, view)

% adds newly detected features to the current map

new_map = old_map;

[V,c] = size(view);

for v=1:V                               % for each feature in a scan
    f = view(v,1);                      % feature id
    if (f ~= 0) & (new_map(f,1) ~= 1)   % this is a new detected feature
        new_map(f,1) = 1;               % detection flag is set
        new_map(f,2:3) = view(v,2:3);
    end
end
