% this function will be used later to see in each iteration at what cost we
% are

function cost = HeatCalCost(X,L,H,tau,params)
    % this quantity is what we aim to minimize
    Z = norm(X-HeatDict(L,tau)*H,"fro").^2;
    fH = params.alpha*sum(sum(abs(H)));
    gL = params.beta*norm(L,"fro").^2;
    cost = Z + fH + gL;
end