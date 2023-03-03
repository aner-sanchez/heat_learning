% tester for different functions

%% HeatGrad.m
% testing to see if we get coherent things
clc, clear, close all

disp("Testing LearnHeat with random parameters")
tau = [2.5;4];
X = randn(3,3);
params = get_params(X,tau);
L = randomLaplacian(params.n);
H = randomSparseH(params,0.2);

[Ldemo,Hdemo,taudemo] = LearnHeat(X,20,0.1,0.01,tau)

%% auxiliar functions

function params = get_params(X,tau)
    params = struct;
    params.n = size(X,1);
    params.m = size(X,2);
    params.S = size(tau,1);
end

function L = randomLaplacian(n)
    L = rand(n);
    L = -(L+L')/2;
    for i=1:n
        L(i,i) = -sum(L(i,:));
    end
end

function H = randomSparseH(params,density)
    H = zeros(params.n*params.S*params.m,1);
    num = round(density*length(H));
    index = randperm(length(H)); index = index(1:num);
    H(index) = 1; H = reshape(H,[params.n*params.S,params.m]);
end