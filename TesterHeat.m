% tester for different functions

%% HeatGrad.m
% testing to see if we get coherent things
clc, clear, close all

disp("Simulation Nº 1")
tau = [2.5;4];
X = randn(3,3);
params = get_params(X,tau);
L = randomLaplacian(params.n);
H = randomSparseH(params,0.2);

[Ldemo,Hdemo,taudemo] = LearnHeat(X,20,0.1,0.01,tau);

%% testing with randgraph_RBF.m and randsignal.m
% currently doesn't learn a laplacian that sums up to 1 by rows
clc, clear

disp("Simulation Nº 2")

% we create a random graph with 10 vertices
L2 = randgraph_RBF(4);

true_tau2 = [2.5;4];
% create 100 random signals
[X2,H2] = randsignal(L2,1,true_tau2);
true_dict2 = HeatDict(L2,true_tau2);

% we will suppose there are two underlying diffusion processes
tau2 = ones(2,1);

% set parameters as paper says
alpha = 1./sqrt(10); beta = 0.1;

[learned_L,learned_H,learned_tau] = LearnHeat(X2,20,alpha,beta,tau2);
X2 - HeatDict(L2,true_tau2)*H2
err2 = testLearnHeat(X2,learned_L,learned_H,learned_tau)
%% testing HeatUpdate.m
clc
L3 = randgraph_RBF(4);
X3 = randsignal(L3,2);
tau3 = ones(2,1);
H3 = rand(8,2);
% X3 - HeatDict(L3,tau3)*H3 % dimensions are ok
params3 = get_params(X3,tau3);
params3.beta = 0.01; params3.alpha = 0.1;
grad_L3 = HeatGradL(X3,L3,H3,tau3,params3);
HeatUpdateL(L3,grad_L3,0.1,params3)

%% auxiliar functions

function error_log = testLearnHeat(X,L,H,tau)
    error_log = X - HeatDict(L,tau)*H;
    error_log(error_log < 1e-5) = 0;
    if (all(all(error_log==0)))
        disp("Perfectly learned!!!!");
    else
        disp('couldnt learn laplacian + signal + tau to perfection');
    end
end

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