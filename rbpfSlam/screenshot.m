%% Probabilistic Machine Learning (CPSC 530A)
%% Instructor: Nando de Freitas
%% Course Project
%% James Cook, Iryna Skrypnyk, Roland Wenzel

% ******************************
% *** Screenshot             ***
% ******************************

function w = screenshot(Mov, frame)

c = clock;
imwrite(frame2im(Mov(frame)), ...
        [date, '-', num2str(c(4)), '.', num2str(c(5)), '.', num2str(fix(c(6))), '.bmp'], ...
        'bmp');
