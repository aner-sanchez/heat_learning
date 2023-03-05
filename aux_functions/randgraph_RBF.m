% creates a random graph by generating n points uniformly distributed in
% [0,1]x[0,1]; it gives the weight W(i,j) = exp(-dist(i,j)^2/2*sigma^2), 

function L = randgraph_RBF(n,sigma2,kappa)
    if nargin < 3; kappa = 0.75; end
    if nargin < 2; sigma2 = 0.5; end
    Z = rand(n,2);
    W = zeros(n);
    for i=1:n
        for j=1:n
            W(i,j) = exp(-(norm(Z(i,:)-Z(j,:)).^2)./(2*sigma2));
        end
    end
    W(W<kappa) = 0;
    L = adj2lap(W);
    L = n.*L/trace(L);
end
%% auxiliar function to convert adjacency matrix to laplacian
function L = adj2lap(W)
    n = size(W,1);
    L = -W;
    for i=1:n
        L(i,i) = sum(W(i,:));
    end
end