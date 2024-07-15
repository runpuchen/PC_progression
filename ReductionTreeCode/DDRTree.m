function [W, Z, stree, Y, history] = DDRTree(X, params)

% X : DxN data matrix
% params.
%       maxIter : maximum iterations
%       eps     : relative objective difference
%       dim     : reduced dimension
%       lambda  : regularization parameter for inverse graph embedding
%       sigma   : bandwidth parameter
%       gamma   : regularization parameter for k-means

[D, N] = size(X);

% initialization
W = pca_projection(X * X', params.dim);
Z = W' * X;

if ~isfield(params,'ncenter')
    K = N; 
    Y = Z(:,1:K);
else
    K = params.ncenter;
    [~, Y] = kmeans(Z',K);
    Y = Y';
end

% main loop
objs = [];
for iter=1:params.maxIter

    
    % Kruskal method to find optimal B
    distsqMU = sqdist(Y,Y);
    

    % stree = graphminspantree(sparse(tril(distsqMU)),'Method','Kruskal');
    % stree = stree + stree';
    G = graph((distsqMU));
    stree = minspantree(G);
    stree = full(adjacency(stree));
    stree = stree+stree';
    B = stree ~= 0;
    L = diag( sum(B,2) ) - B;
    
    % compute R using mean-shift update rule
    distZY = sqdist(Z,Y);
    min_dist = repmat(min(distZY,[],2),1,K);
    tmp_distZY = distZY - min_dist;
    tmp_R = exp(-tmp_distZY ./ params.sigma);
    R = tmp_R ./ repmat( sum(tmp_R,2), 1, K );
    Gamma = diag( sum(R) );
    
    % termination condition
    obj1 = - params.sigma * sum( log( sum( exp(- tmp_distZY./params.sigma) ,2) ) ...
    - min_dist(:,1) ./ params.sigma );
    objs(iter) = (norm(X-W*Z))^2 + params.lambda .* trace( Y * L * Y' )...
        + params.gamma * obj1;
    fprintf('iter=%d obj = %f\n',iter,objs(iter));
    
    history.W{iter} = W;
    history.Z{iter} = Z;
    history.Y{iter} = Y;
    history.stree{iter} = stree;
    history.R{iter} = R;
    
    if iter >1 
        if abs(objs(iter) - objs(iter-1))/abs(objs(iter-1)) < params.eps
            break;
        end
    end 
    
    % compute low dimension projection matrix
    tmp = R / ( ((params.gamma +1) / params.gamma) .* ( (params.lambda / params.gamma)  .* L + Gamma) - R' * R);
    Q = 1/(params.gamma +1) .* ( eye(N,N) + tmp * R' );
    C = X * Q;
    tmp1 = C * X';
    W = pca_projection( (tmp1 + tmp1')./2, params.dim);
    Z = W' * C;
    Y = Z * R / (params.lambda / params.gamma .* L + Gamma);
end

history.objs = objs;