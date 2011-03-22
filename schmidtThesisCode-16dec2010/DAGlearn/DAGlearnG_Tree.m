function [adj] = DAGlearn2_ChowLiu(X,scoreType,A)
% Optimal tree for Gaussian CPDs based on BIC or Validation

if nargin < 3
    A = [];
end

[nSamples,nNodes] = size(X);

% Compute unary scores
options.Display = 0;
for n1 = 1:nNodes
    %fprintf('Processing node %d...\n',n);
    if isempty(A)
        intInd = [];
    else
        intInd = A(:,n1)~=0;
    end
    
    if scoreType == 0
        if isempty(intInd)
            ysub = X(:,n1);
        else
            nSamplesSub = sum(~intInd);
            ysub = X(~intInd,n1);
        end
        n = size(ysub,1);
        sigma2 = sum((ysub).^2)/n;
        nll = n*log(sqrt(sigma2)) + (n/2)*log(2*pi) + (norm(ysub)^2)/(2*sigma2);
        score1(n1,1) = 2*nll;
    else
        trainNdx = [1:nSamples]' <= ceil(nSamples/2);
        if isempty(intInd)
            ysub = X(trainNdx,n1);
        else
            ysub = X(trainNdx & ~intInd,n1);
        end
        n = size(ysub,1);
        sigma2 = sum((ysub).^2)/n;
        if isempty(intInd)
            n = sum(~trainNdx);
            ysub = X(~trainNdx,n1);
        else
            n = sum(~trainNdx & ~intInd);
            ysub = X(~trainNdx & ~intInd,n1);
        end
        n = length(ysub);
        score1(n1,1) = n*log(sqrt(sigma2)) + (n/2)*log(2*pi) + (norm(ysub)^2)/(2*sigma2);
    end
end

% Compute conditional scores
for n1 = 1:nNodes
    fprintf('Processing node %d...\n',n1);
    for n2 = 1:nNodes
        if n1 == n2
            continue
        end
        %fprintf('Processing nodes (%d | %d)...\n',n1,n2);
        
        if isempty(A)
            intInd = [];
        else
            intInd = A(:,n1)~=0;
        end
        
        if scoreType == 0
            if isempty(intInd)
                Xsub = X(:,n2);
                ysub = X(:,n1);
            else
                Xsub = X(~intInd,n2);
                ysub = X(~intInd,n1);
            end
            w = Xsub\ysub;
            n = size(Xsub,1);
            sigma2 = sum((Xsub*w - ysub).^2)/n;
            nll = n*log(sqrt(sigma2)) + (n/2)*log(2*pi) + (norm(Xsub*w-ysub)^2)/(2*sigma2);
            score2(n1,n2) = 2*nll + length(w)*log(nSamples);
        else
           if isempty(intInd)
               Xsub = X(trainNdx,n2);
               ysub = X(trainNdx,n1);
           else
               Xsub = X(trainNdx & ~intInd,n2);
               ysub = X(trainNdx & ~intInd,n1);
           end
           w = Xsub\ysub;
            n = size(Xsub,1);
            sigma2 = sum((Xsub*w - ysub).^2)/n;
           if isempty(intInd)
               Xsub = X(~trainNdx,n2);
               ysub = X(~trainNdx,n1);
           else
               Xsub = X(~trainNdx & ~intInd,n2);
               ysub = X(~trainNdx & ~intInd,n1);
           end
           n = size(Xsub,1);
           score2(n1,n2) = n*log(sqrt(sigma2)) + (n/2)*log(2*pi) + (norm(Xsub*w-ysub)^2)/(2*sigma2);
        end
    end
end

if 0 % Show weight matrix
    for n1 = 1:nNodes
        for n2 = 1:nNodes
            weights(n1,n2) = score1(n1) - score2(n1,n2);
        end
    end
    weights
    max(max(abs(weights-weights')))
    pause
end

if isempty(A) || sum(A(:)) == 0
    % In observational case, have score-equivalence so just use min-weight
    % spanning tree (currently, this gives correct adjacency but score
    % calculation is incorrect)
    
    %% Weights are conditionals-unconditionals
    edgeEnds = zeros(0,3);
    for n1 = 1:nNodes
        for n2 = n1+1:nNodes
            edgeEnds(end+1,:) = [n1 n2 score2(n1,n2) - score1(n1)];
        end
    end
    
    %% Solve
    E = minSpan(nNodes,edgeEnds);
    
    %% Make set of selected edges
    edges = zeros(sum(E),2);
    e2 = 1;
    for e = 1:length(E)
        if E(e)==1
            edges(e2,:) = edgeEnds(e,1:2);
            e2 = e2+1;
        end
    end
    
    %% Make adjacency matrix
    adj = zeros(nNodes);
    for e = 1:length(E)
        if E(e)==1
            adj(edgeEnds(e,1),edgeEnds(e,2)) = 1;
        end
    end
    %drawGraph(adj)
    
    %% Compute score
    if 0
        score = 0;
        for n = 1:nNodes
            parent = find(adj(:,n))
            if isempty(parent)
                score = score + score1(n);
            else
                score = score + score2(n,parent);
            end
        end
    end
else
    % With interventional data, need to compute maximal branching
    
    %% Weights are conditionals-unconditionals
    edgeEnds = zeros(0,3);
    for n1 = 1:nNodes
        for n2 = 1:nNodes
            if n1 == n2
                continue
            end
            edgeEnds(end+1,:) = [n2 n1 score1(n1) - score2(n1,n2)];
        end
    end
    
    %% Solve
    fprintf('Running Edmonds algorithm...\n');
    GT = edmonds(1:nNodes,edgeEnds);
    selected=reconstruct_2(GT,0);
    edgeEnds = edgeEnds(selected,1:2);
    
    %% Make adjacency matrix
    adj = zeros(nNodes);
    for e = 1:size(edgeEnds,1)
        adj(edgeEnds(e,1),edgeEnds(e,2)) = 1;
    end
    %drawGraph(adj)
    
    %% Compute score
    if 0
        score = 0;
        for n = 1:nNodes
            parent = find(adj(:,n))
            if isempty(parent)
                score = score + score1(n);
            else
                score = score + score2(n,parent);
            end
        end
    end
end