clear all;

randn('seed',0);
rand('seed',0);
fprintf('Generating data...\n');
% Data are generated
I=double(imread('data/lena.png'))/255;
% extract 10 x 10 patches
X=im2col(I,[10 10],'sliding');
X=X(:,1:10:end);
per=randperm(size(X,2));
D=X(:,per(1:200));
X=X-repmat(mean(X),[size(X,1) 1]);
X=X ./ repmat(sqrt(sum(X.^2)),[size(X,1) 1]);
 D=D-repmat(mean(D),[size(D,1) 1]);
 D=D ./ repmat(sqrt(sum(D.^2)),[size(D,1) 1]);
 X=X(:,1:1000);

% parameter of the optimization procedure are chosen
param.lambda1=0.15;
param.lambda=param.lambda1;
param.lambda2=0.05;
param.lambda3=0.15;
param.numThreads=4; % number of processors/cores to use; the default choice is -1
param.mode=2;       % penalized formulation
param.tol=0.001;
param.itermax1=50;
param.accelerated=false;
param.rho=0.1;


fprintf('Processing data...\n');
tic
alpha3=mexNesterov(X,D,zeros(size(D,2),size(X,2)),param);
t=toc;
toc

fprintf('%f signals processed per second for Nesterov non accelerated \n',size(X,2)/t);
E=mean(0.5*sum((X-D*alpha3).^2)+param.lambda1*sum(abs(alpha3))+0.5*param.lambda2*sum(alpha3.^2));
fprintf('Objective function for Nesterov 1: %g\n',E);



param.accelerated=true;

tic
alpha3=mexNesterov(X,D,zeros(size(alpha3)),param);
t=toc;
toc

fprintf('%f signals processed per second for Nesterov accelerated \n',size(X,2)/t);
E=mean(0.5*sum((X-D*alpha3).^2)+param.lambda1*sum(abs(alpha3))+0.5*param.lambda2*sum(alpha3.^2));
fprintf('Objective function for Nesterov 1: %g\n',E);

