function [W, Z, stree, history] = DRTree(X, params)

% X : DxN data matrix
% params.
%       maxIter : maximum iterations
%       eps     : relative objective difference
%       dim     : reduced dimension
%       lambda  : regularization parameter for inverse graph embedding

[D, N] = size(X);

% initialization
W = pca_projection(X * X', params.dim);
Z = W' * X;

% main loop
objs = [];
for iter=1:params.maxIter
    
    % Kruskal method to find optimal B
    distsqMU = sqdist(Z,Z);
    stree = graphminspantree(sparse(tril(distsqMU)),'Method','Kruskal');
    stree = stree + stree';
    B = stree ~= 0;
    L = diag( sum(B,2) ) - B;
    
    % termination condition
    objs(iter) = (norm(X-W*Z))^2 + params.lambda * trace( Z * L * Z' );
    fprintf('iter=%d obj = %f\n',iter,objs(iter));
    
    history.W{iter} = W;
    history.Z{iter} = Z;
    history.stree{iter} = stree;
    
    if iter >1 
        if abs(objs(iter) - objs(iter-1))/abs(objs(iter-1)) < params.eps
            break;
        end
    end 
    
    % compute W and Z
    tmp = X / (eye(N,N) + params.lambda .* L);
    C = tmp * X';
    W = pca_projection( (C+C')./2, params.dim);
    Z = W' * tmp;  
end

history.objs = objs;
