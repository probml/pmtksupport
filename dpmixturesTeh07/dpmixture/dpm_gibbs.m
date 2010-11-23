function dpm = dpm_gibbs(dpm,numiter);
% run numiter number of iterations of gibbs sampling in the DP mixture

KK = dpm.KK; % number of active clusters
NN = dpm.NN; % number of data items
aa = dpm.aa; % alpha parameter
qq = dpm.qq; % row cell vector of mixture components
xx = dpm.xx; % row cell vector of data items
zz = dpm.zz; % row vector of cluster indicator variables
nn = dpm.nn; % row vector of number of data items per cluster

for iter = 1:numiter
  % in each iteration, remove each data item from model, then add it back in
  % according to the conditional probabilities.

  for ii = 1:NN % iterate over data items ii

    % remove data item xx{ii} from component qq{kk}
    kk = zz(ii); % kk is current component that data item ii belongs to
    nn(kk) = nn(kk) - 1; % subtract from number of data items in component kk
    qq{kk} = delitem(qq{kk},xx{ii}); % subtract data item sufficient statistics

    % delete active component if it has become empty
    if nn(kk) == 0, 
      %fprintf(1,'del component %3d. K=%3d\n',find(nn==0),KK-sum(nn==0));
      KK = KK - 1;
      qq(kk) = [];
      nn(kk) = [];
      idx = find(zz>kk);
      zz(idx) = zz(idx) - 1;
    end

    % compute conditional probabilities pp(kk) of data item ii
    % belonging to each component kk
    % compute probabilities in log domain, then exponential
    pp = log([nn aa]);
    for kk = 1:KK+1
      pp(kk) = pp(kk) + logpredictive(qq{kk},xx{ii});
    end
    pp = exp(pp - max(pp)); % -max(p) for numerical stability
    pp = pp / sum(pp);

    % choose component kk by sampling from conditional probabitilies
    uu = rand;
    kk = 1+sum(uu>cumsum(pp));

    % instantiates a new active component if needed
    if kk == KK+1
      %fprintf(1,'add component %3d. K=%3d\n',kk,KK+1);
      KK = KK + 1;
      nn(kk) = 0;
      qq(kk+1) = qq(kk);
    end

    % add data item xx{ii} back into model (component qq{kk})
    zz(1,ii) = kk; 
    nn(1,kk) = nn(1,kk) + 1; % increment number of data items in component kk
    qq{1,kk} = additem(qq{1,kk},xx{ii}); % add sufficient stats of data item

  end
end

% save variables into dpm struct
dpm.qq = qq;
dpm.zz = zz;
dpm.nn = nn;
dpm.KK = KK;
