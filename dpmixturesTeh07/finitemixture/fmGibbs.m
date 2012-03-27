function fm = fm_gibbs(fm,numiter);
% run numiter number of iterations of gibbs sampling in the finite mixture

KK = fm.KK;
NN = fm.NN;
aa = fm.aa;
qq = fm.qq;
xx = fm.xx;
zz = fm.zz;
nn = fm.nn;

for iter = 1:numiter
  % in each iteration, remove each data item from model, then add it back in.

  for ii = 1:NN

    % remove data item xx{ii} from component qq{kk}
    kk = zz(ii);
    nn(kk) = nn(kk) - 1;
    qq{kk} = delitem(qq{kk},xx{ii});

    % compute probabilities pp(kk) of each component kk
    pp = log(aa/KK + nn);
    for kk = 1:KK
      pp(kk) = pp(kk) + logpredictive(qq{kk},xx{ii});
    end
    pp = exp(pp - max(pp));
    pp = pp / sum(pp);

    % choose component kk
    uu = rand;
    kk = 1+sum(uu>cumsum(pp));

    % add data item xx{ii} back into model (component qq{kk})
    zz(ii) = kk;
    nn(kk) = nn(kk) + 1;
    qq{kk} = additem(qq{kk},xx{ii});

  end
end

fm.qq = qq;
fm.zz = zz;
fm.nn = nn;

