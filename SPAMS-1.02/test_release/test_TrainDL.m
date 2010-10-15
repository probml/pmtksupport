clear all; 
%param.approx=4.0;  EXPERIMENTAL OPTION

I=double(imread('data/lena.png'))/255;
% extract 8 x 8 patches
X=im2col(I,[8 8],'sliding');
X=X-repmat(mean(X),[size(X,1) 1]);
X=X ./ repmat(sqrt(sum(X.^2)),[size(X,1) 1]);


param.K=100;  % learns a dictionary with 100 elements
param.lambda=0.15;
param.numThreads=4; % number of threads
param.batchsize=512;

param.iter=100;  % let us see what happens after 100 iterations.


%%%%%%%%%% FIRST EXPERIMENT %%%%%%%%%%%
tic
D = mexTrainDL(X,param);
t=toc;
fprintf('time of computation for Dictionary Learning: %f\n',t);

param.approx=0;
fprintf('Evaluating cost function...\n');
alpha=mexLasso(X,D,param);
R=mean(0.5*sum((X-D*alpha).^2)+param.lambda*sum(abs(alpha)));
ImD=displayPatches(D);
subplot(1,3,1);
imagesc(ImD); colormap('gray');
fprintf('objective function: %f\n',R);


return;

fprintf('*********** SECOND EXPERIMENT ***********\n');

%%%%%%%%%% SECOND EXPERIMENT %%%%%%%%%%%
% Train on half of the training set, then retrain on the second part
X1=X(:,1:floor(size(X,2)/2));
X2=X(:,floor(size(X,2)/2):end);
param.iter=50;
tic
[D model] = mexTrainDL(X1,param);
t=toc;
fprintf('time of computation for Dictionary Learning: %f\n',t);
fprintf('Evaluating cost function...\n');
alpha=mexLasso(X,D,param);
R=mean(0.5*sum((X-D*alpha).^2)+param.lambda*sum(abs(alpha)));
fprintf('objective function: %f\n',R);
tic
% Then reuse the learned model to retrain a few iterations more.
param2=param;
param2.D=D;
[D model] = mexTrainDL(X2,param2,model);
%[D] = mexTrainDL(X,param);
t=toc;
fprintf('time of computation for Dictionary Learning: %f\n',t);
fprintf('Evaluating cost function...\n');
alpha=mexLasso(X,D,param);
R=mean(0.5*sum((X-D*alpha).^2)+param.lambda*sum(abs(alpha)));
fprintf('objective function: %f\n',R);

%%%%%%%% FOURTH EXPERIMENT %%%%%%%%%%%%%%%%%%%%%%
% let us add sparsity to the dictionary itself

fprintf('*********** FOURTH EXPERIMENT ***********\n');
param.modeParam=0;
param.iter=100;
param.gamma1=0.3;
param.modeD=1;
tic
[D] = mexTrainDL(X,param);
t=toc;
fprintf('time of computation for Dictionary Learning: %f\n',t);
fprintf('Evaluating cost function...\n');
alpha=mexLasso(X,D,param);
R=mean(0.5*sum((X-D*alpha).^2)+param.lambda*sum(abs(alpha)));
fprintf('objective function: %f\n',R);
tic
subplot(1,3,2);
ImD=displayPatches(D);
imagesc(ImD); colormap('gray');


fprintf('*********** FIFTH EXPERIMENT ***********\n');
param.modeParam=0;
param.iter=100;
param.gamma1=0.3;
param.gamma2=0.3;
param.modeD=3;
tic
[D] = mexTrainDL(X,param);
t=toc;
fprintf('time of computation for Dictionary Learning: %f\n',t);
fprintf('Evaluating cost function...\n');
alpha=mexLasso(X,D,param);
R=mean(0.5*sum((X-D*alpha).^2)+param.lambda*sum(abs(alpha)));
fprintf('objective function: %f\n',R);
tic
subplot(1,3,3);
ImD=displayPatches(D);
imagesc(ImD); colormap('gray');




