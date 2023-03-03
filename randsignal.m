% generates m graph signals by linearly combining three random atoms from
% the dictionary with random coefficients drawn from a Gaussian
% distribution with zero mean and unit variance

function [X,H] = randsignal(L,m,tau)
    n = size(L,1);
    if nargin < 3; tau = [2.5;4]; end
    s = size(tau,1);
    D = HeatDict(L,tau);
    X = []; H = [];
    for i=1:m
        randindex = randperm(size(D,2),3);
        randcoeff = randn(3,1);
        h = zeros(n*s,1); h(randindex) = randcoeff;
        x = D(:,randindex)*randcoeff;
        X = [X,x]; H = [H,h];
    end
end