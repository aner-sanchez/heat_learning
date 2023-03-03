% Creates dictionary from a vector tau [tau1,...,tauS] and a laplacian
% matrix L, D = [expm(-tau1*L),...,expm(-tauS*L)]; see Algorithm 1
% https://arxiv.org/abs/1611.01456

function D = HeatDict(L,tau,D)
    assert(size(tau,2)==1,['vector tau of diffusion processes not provided...' ...
        'in column form to HeatDict algorithm'])
    assert(size(L,1)==size(L,2),'Laplacian is not squared');
    S = size(tau,1); N = size(L,1);
    if nargin < 3; D = zeros(N,N*S);
    else; assert(size(D,2)/size(D,1) == S,['dictionary size provided is not in accordance...' ...
            'to tau and L provided']);
    end
    for i=1:S
        D(:,1+(i-1)*N:i*N) = expm(-tau(i,1).*L);
    end
end