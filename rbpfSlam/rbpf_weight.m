%% Probabilistic Machine Learning (CPSC 530A)
%% Instructor: Nando de Freitas
%% Course Project
%% James Cook, Iryna Skrypnyk, Roland Wenzel

% ******************************
% *** RBPF_Weight            ***
% ******************************

function w = rbpf_weight(curr_view, proposed_view, sigma)

w = 1;
[F,c] = size(curr_view);

% compute joint probability (assume independency)
for f = 1:F     % for each feature
    if curr_view(f,1) ~= 0 & curr_view(f,2) ~= 0
        w = w * multigaussian(curr_view(f,:)', proposed_view(f,:)', sigma(f,f)*eye(2));
    end
end

% ******************************
% *** Multigaussian          ***
% ******************************

function y = multigaussian(x, mu, sigma)

% Multigaussian distribution (up to a normalization constant)
y = exp(-0.5 * (x - mu)' * inv(sigma) * (x - mu)) + 1e-99; % (det(sigma)^(-0.5))*...
