%%%%%%%%%%%%%%%%% State Definition %%%%%%%%%%%%%
%
%    S_k = A_k S_{k-1} + N(0 , R_k)
%
%    S_k = (x_k , x'_k , y_k , y'_k , H_k^x , H_k^y , theta_k)
%
%
%          |1   delta_k       0        0       0      0     0|
%          |0         1       0        0       0      0     0|
%          |0         0       1  delta_k       0      0     0|
% A  =     |0         0       0        1       0      0     0|
%          |0         0       0        0       1      0     0|
%          |0         0       0        0       0      1     0|
%          |0         0       0        0       0      0     1|
%
% R_k  =   |R_y       0   |
%          |0         R_e |
%
%                 |delta_k^3/3   delta_k^2/2       0        0 |
% R_y  = sigma_y  |delta_k^2/2       delta_k       0        0 |, % Kinematic of
%                 |0         0       delta_k^3/3   delta_k^2/2|  % the ellipsoid
%                 |0         0       delta_k^2/2       delta_k|

%          |sigma_Hx^2        0          0|
% R_e  =   |0        sigma_Hy^2          0|, % Kinematic of
%          |0         0      sigma_theta^2|  % the ellipsoid




clear all

close all


video             = 'Second.avi';

info              = aviinfo(video);

offset_frame      = 80;%80

%nb_frame          = info.NumFrames - offset_frame -10;

nb_frame          = 400;

dim_x             = info.Width;

dim_y             = info.Height;



N                 = 250;        % Number of particules

N_threshold       = 6.*N/10;    % Redistribution threshold

delta             = 0.7;

%%%%% Color Cue parameters %%%%%%%%

Npdf              = 120;       % Number of samples to draw inside ellipse to evaluate color histogram

Nx                = 6;         % Number of bins in first color dimension (R or H)

Ny                = 6;         % Number of bins in second color dimension (G or S)

Nz                = 6;         % Number of bins in third color dimension (B or V)

sigma_color       = 0.20;      % Measurement Color noise

nb_hist           = 256;




range             = 1;

pos_index         = [1 , 3];

ellipse_index     = [5 , 6 , 7];

d                 = 7;

M                 = Nx*Ny*Nz;

vect_col          = (0:range/(nb_hist - 1):range);



%%%%%% Target Localization for computing the target distribution %%%%

yq                = [185 ; 100];

eq                = [14 ; 20 ; pi/3];

%%%%%% Initialization distribution initialization %%%%

Sk                = zeros(d , 1);

Sk(pos_index)     = yq;

Sk(ellipse_index) = eq;


% Initial State covariance %


sigmax1           = 60;      % pixel %

sigmavx1          = 1;      % pixel / frame %

sigmay1           = 60;      % pixel  %

sigmavy1          = 1;      % pixel / frame %


sigmaHx1          = 4;    % pixel %

sigmaHy1          = 4;    % pixel %

sigmatheta1       = 8*(pi/180); % rad/frame %


% State Covariance %

% a) Position covariance %

sigmay            = sqrt(0.35);


% b) ellipse covariance %

sigmaHx           = 0.1;                % pixel %

sigmaHy           = 0.1;                % pixel %

sigmatheta        = 3.0*(pi/180);       % rad/frame %


%%%%%%%%%%%%%%%%%%%% State transition matrix %%%%%%%%%%%%%%%%%%%%%%

A                 = [1 delta 0 0 0 0 0 ; 0 1 0 0 0 0 0 ; 0 0 1 delta 0 0 0; 0 0 0 1 0 0 0 ; 0 0 0 0 1 0 0 ; 0 0 0 0 0 1 0 ; 0 0 0 0 0 0 1];

By                = [1 0 0 0 0 0 0 ; 0 0 1 0 0 0 0];

Be                = [0 0 0 0 1 0 0 ; 0 0 0 0 0 1 1 ; 0 0 0 0 0 0 1 ];

%%%%%% Initial State Covariance %%%%%

R1                = diag([sigmax1 , sigmavx1 , sigmay1 , sigmavy1 , sigmaHx1 , sigmaHy1 , sigmatheta1].^2);

%%%%%% State Covariance %%%%%

Rk                = zeros(d , d);

Ry                = (sigmay^2)*[delta^3/3 delta^2/2 0 0 ; delta^2/2 delta 0 0 ; 0 0 delta^3/3 delta^2/2 ; 0 0 delta^2/2 delta];

Re                = [sigmaHx.^2 0 0 ; 0 sigmaHy.^2 0 ; 0 0 sigmatheta.^2];


Rk(1 : 4  , 1 : 4)= Ry;

Rk(5 : d , 5 : d) = Re;

Ck                = chol(Rk)';


%%%%%%%% Memory Allocation %%%%%%%

ON                = ones(1 , N);

Od                = ones(d , 1);

Smean             = zeros(d , nb_frame);

Pcov              = zeros(d , d , nb_frame);

N_eff             = zeros(1 , nb_frame);

cte               = 1/N;

cteN              = cte(1 , ON);

w                 = cteN;


compteur          = 0;

cte1_color        = 1/(2*sigma_color*sigma_color);

cte2_color        = (1/(sqrt(2*pi)*sigma_color));



%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Target Distribution %%%%%%%%%%%%%%

warning('off','MATLAB:aviread:FunctionToBeRemoved')

mov               = aviread(video , offset_frame + 1);

im                = mov.cdata;

Z                 = double(im);

im                = rgb2hsv_mex(Z);



C1                = cumsum(histc(reshape(im(: , : , 1) , dim_x*dim_y , 1) , vect_col))/(dim_x*dim_y);

C2                = cumsum(histc(reshape(im(: , : , 2) , dim_x*dim_y , 1) , vect_col))/(dim_x*dim_y);

C3                = cumsum(histc(reshape(im(: , : , 3) , dim_x*dim_y , 1) , vect_col))/(dim_x*dim_y);


i1                = sum(C1(: , ones(1 , Nx)) < repmat((0:1/(Nx - 1) : 1) , nb_hist  , 1));

i2                = sum(C2(: , ones(1 , Ny)) < repmat((0:1/(Ny - 1) : 1) , nb_hist  , 1));

i3                = sum(C3(: , ones(1 , Nz)) < repmat((0:1/(Nz - 1) : 1) , nb_hist  , 1));


edge1             = [0 , vect_col(i1(2 : end)) , range];

edge2             = [0 , vect_col(i2(2 : end)) , range];

edge3             = [0 , vect_col(i3(2 : end)) , range];

q                 = pdfcolor_ellipserand(im , yq , eq , Npdf , edge1 , edge2 , edge3);

Q                 = q(: , ON);

fig1              = figure(1);

image(mov.cdata);

set(gca , 'drawmode' , 'fast');

set(gcf , 'doublebuffer','on');

set(gcf , 'renderer' , 'zbuffer');



%%%%%%%%%%%% Particle initialisation %%%%%%%%

Sk                = Sk(: , ON) + chol(R1)'*randn(d , N);

%%%%%%%%%%%% Main Loop %%%%%%%%%%%%%%%%%%%%%%

for k = 1 : nb_frame ;


    disp(sprintf('Frames = %d/%d' , k , nb_frame ));

    Sk               = A*Sk + Ck*randn(d , N);

    mov              = aviread(video , offset_frame + k);

    im               = (mov.cdata);

    Z                = double(im);

    im               = rgb2hsv_mex(Z);


    yk               = Sk(pos_index , :);     %yk                = By*Sk;

    ek               = Sk(ellipse_index , :); %ek                = Be*Sk;

    %%%%%%%%%%  Color Likelihood %%%%%%%%%

    [py , zi , yi]   = pdfcolor_ellipserand(im , yk , ek , Npdf , edge1 , edge2 , edge3);

    rho_py_q         = sum(sqrt(py.*Q));

    likelihood_color = cte2_color*exp((rho_py_q - 1)*cte1_color);

    w                = w.*likelihood_color;

    w                = w/sum(w);

    %--------------------------- 6) MMSE estimate & covariance -------------------------

    [Smean(: , k) , Pcov(: , : , k)] = part_moment(Sk , w);

    %--------------------------- 7) Particles redistribution ? if N_eff < N_threshold -------------------------

    N_eff(k)                         = 1./sum(w.*w);

    if (N_eff(k) < N_threshold)

        compteur              = compteur + 1;

        indice_resampling     = particle_resampling(w);

        % Recopie des particules selon le tirage des indices précédents

        Sk                    = Sk(: , indice_resampling);

        w                     = cteN;

    end


    %%%%%%%%%%%%% Display %%%%%%%%%%%%%%%


    fig1              = figure(1);

    image(mov.cdata);

    title(sprintf('N = %6.3f/%6.3f, Frame = %d, Redistribution =%d' , N_eff(k) , N_threshold , k , compteur))


    ind_k             = (1 : k);

    hold on

    ykmean            = Smean(pos_index , k);

    ekmean            = Smean(ellipse_index, k);

    [xmean , ymean]   = ellipse(ykmean , ekmean);


    plot(xmean , ymean , 'g' , 'linewidth' , 3)

    plot(Smean(pos_index(1) , ind_k) , Smean(pos_index(2) , ind_k) , 'r' , 'linewidth' , 2)

    plot(Sk(pos_index(1) , :) , Sk(pos_index(2) , :) , 'b+');

    hold off

    pause
end





figure(3)

plot(Smean(pos_index(1) , :) , Smean(pos_index(2) , :) , 'linewidth' , 2)

axis([1 , dim_x , 1 , dim_y ]);

axis ij

grid on

title('Helicopter trajectory');


figure(4)

plot((1 : nb_frame) , Smean(ellipse_index(end) , :) , 'linewidth' , 2)

axis([1 , nb_frame , -pi , pi ]);

xlabel('Frames k');

ylabel('\theta');

grid on

title('Ellipse angle versus frames');



figure(5);

h = slice(reshape(q , Nx , Ny , Nz) , (1 : Nx) , (1 : Ny) , (1 : Nz));

colormap(flipud(cool));

alpha(h , 0.1);

brighten(+0.5);

title('3D Histogram of the Target distribution');

xlabel('Bin x');

ylabel('Bin y');

zlabel('Bin z');

colorbar

cameramenu;


figure(6);

hold on


plot(vect_col , C1 , vect_col , C2, vect_col , C3);

plot(edge1 , C1([1 , i1(2 : end) , nb_hist]) , '+' , edge2 , C2([1 , i2(2 : end) , nb_hist])   ,'*' , edge3 , C3([1 , i3(2 : end) , nb_hist]) ,'p')

hold off

ylabel('HSV CDF');



