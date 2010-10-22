function [G,local,lambda,exp_local,exp_lambda,xs,s_pair,Fields] = tst_potts(graph,nVals,nPerSide,inter_scheme,inter_strength,field_strength)
%nPerSide =5;
%inter_scheme =3;
%inter_strength = 1;
%field_strength = 0.05;
%nVals =3;

%graph = 'grid';
%graph = 'cycle';
%graph = 'clique';
switch graph
  case 'grid2'
    [G,xs,s_pair] =  make_grid_graph(nPerSide,1,inter_scheme,inter_strength,field_strength);
  case 'grid'
    [G,xs,s_pair] =  make_grid_graph(nPerSide,nPerSide,inter_scheme,inter_strength,field_strength);
  case 'tree'
    % Make a tree
    G.Edges = zeros(5);
    G.Edges(1,[2 3]) = 1;
    G.Edges(2,[4 5])=1;
    G.Edges = G.Edges+G.Edges';
    G.Weights = G.Edges.*rand(5);
    G.Fields = rand(5,1);
  case 'cycle'
    % Make a cycle
    G.Edges = diag(ones(nPerSide-1,1),1);
    G.Edges(1,end)=1;
    G.Edges = G.Edges+G.Edges';
    G.Fields = (rand(nPerSide,1)-0.5)*field_strength;
    G.Fields(3)=0.1;
    switch inter_scheme
        case 1
            ww = (rand(size(G.Edges)))*inter_strength;
        case 2
            ww =(-rand(size(G.Edges)))*inter_strength;
        case 3
            ww = 2.*(rand(size(G.Edges))-0.5)*inter_strength;
    end
    G.Weights = G.Edges.*ww;
    G.Weights = (G.Weights + G.Weights')/2;
    G.Weights(1,3)=eps;
    G.Weights(2,3)=eps;
    G.Weights(3,1)=eps;
    G.Weights(3,2)=eps;
  case 'clique'
    G.Edges = ones(nPerSide)-eye(nPerSide);
    G.Fields = (rand(nPerSide,1)-0.5)*field_strength;
    switch inter_scheme
        case 1
            ww = (rand(size(G.Edges)))*inter_strength;
        case 2
            ww =(-rand(size(G.Edges)))*inter_strength;
        case 3
            ww = 2.*(rand(size(G.Edges))-0.5)*inter_strength;
    end
    G.Weights = G.Edges.*ww;
    G.Weights = (G.Weights + G.Weights')/2;
    G.Fields = G.Fields*0+1;
%    w = [0 1 2 3; 3 0 2 1; 1 2 0 3; 2 2 1 0];
    w = round(rand(nPerSide)*20); w = w-diag(diag(w));
    G.Weights = (w+w').*sign(ww);
    G.ArrWeights=[];
end


nVertices = size(G.Edges,1);
degs = sum(G.Edges,2);
lambda=[]; local=[]; exp_lambda=[]; exp_local = [];
mn = Inf;
mx = -Inf;
%if ~strcmp(graph,'cycle')
Fields = 2*(rand(nVertices,nVals)-0.5)*field_strength;
%else
%  Fields = G.Fields;
%xsend

%Fields =(ones(nVertices,nVals));
%Fields = G.Fields;
for v1=1:nVertices
    local{v1} = zeros(nVals,1);
    exp_local{v1} = exp(local{v1});
    for v2=1:nVertices
        if v1==v2 | ~G.Edges(v1,v2)
            continue;
        end
        w = G.Weights(v1,v2);
        if nVals==2 & strcmp(graph,'cycle')
            f1 = G.Fields(v1)/degs(v1);
            f2 = G.Fields(v2)/degs(v2);
            lambda{v1,v2} = [w-f1-f2 -w-f1+f2; -w+f1-f2 w+f1+f2]/2;
	    exp_lambda{v1,v2} = exp(lambda{v1,v2});;
            mn = min(mn,min(lambda{v1,v2}(:)));
            mx = max(mx,max(lambda{v1,v2}(:)));	    
        else
            for i1=1:nVals
                for i2=1:nVals
                    f1 = Fields(v1,i1)/degs(v1);
                    f2 = Fields(v2,i2)/degs(v2);
                    lambda{v1,v2}(i1,i2) = w*(i1==i2)+f1+f2;
                    exp_lambda{v1,v2} = exp(lambda{v1,v2});
                end
            end
            mn = min(mn,min(lambda{v1,v2}(:)));
            mx = max(mx,max(lambda{v1,v2}(:)));	    
        end
    end
end

for v1=1:nVertices
    for v2=v1+1:nVertices
        w = G.Weights(v1,v2);
        f1 = G.Fields(v1);
        f2 = G.Fields(v2);
%        lambda{v1,v2} = lambda{v1,v2}-mn;
        lambda{v1,v2} = lambda{v2,v1}';
	exp_lambda{v1,v2} = exp_lambda{v2,v1}';
    end
end

%map_msgs = mymap(G.Edges,lambda);
return
msgs = map_msgs;
[row_inds,col_inds] = find(triu(G.Edges));
for i=1:length(row_inds)
    ri = row_inds(i);
    ci = col_inds(i);
    t = lambda{ri,ci}-repmat(msgs{ri,ci},2,1)-repmat(msgs{ci,ri},1,2);
    t
    pause
end


