% tester for different functions

%% HeatGrad.m
% testing to see if we get coherent things
clc, clear, close all

disp("Simulation Nº 1")
tau = [2.5;4];
X = randn(10,3);
params = get_params(X,tau);
L = randomLaplacian(params.n);
H = randomSparseH(params,0.2);

[Ldemo,Hdemo,taudemo] = LearnHeat(X,100,0.1,0.01,tau);

%% testing with randgraph_RBF.m and randsignal.m
% currently doesn't learn a laplacian that sums up to 1 by rows
clc, clear

disp("Simulation Nº 2")
disp("..............................")

n = 5; m = 2;

% we create a random graph with 10 vertices
L2 = randgraph_RBF(n);

true_tau2 = [2.5;4];
% create 100 random signals
[X2,H2] = randsignal(L2,m,true_tau2);
true_dict2 = HeatDict(L2,true_tau2);

% we will suppose there are two underlying diffusion processes
tau2 = ones(2,1);

% set parameters as paper says
alpha = 1./sqrt(10); beta = 0.1;

[learned_L,learned_H,learned_tau] = LearnHeat(X2,20,alpha,beta,tau2,0);
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

%% testing backtrackingL
clc,clear
disp("Testing backtrackingL")
n = 5; m = 20;
L2 = randgraph_RBF(n);
true_tau2 = [2.5;4];
[X2,H2] = randsignal(L2,m,true_tau2);
true_dict2 = HeatDict(L2,true_tau2);
alpha = 1./sqrt(10); beta = 0.1;
params = struct("beta",beta,"n",n,"S",2); lip = struct("d",1.1); gamma = struct("d",1.1);
L_new = backtrackingL(X2,L2,H2,true_tau2,params,lip,gamma);
%% testing python program
clc, clear
L = [1 0 -1 0;
    0 2 -1 -1;
    -1 -1 2 0;
    0 -1 0 1];
tau = [2.5,4];
H = [ 0.       0.03763  0.       0.3431 ;
  0.       0.       0.20246  0.     ;
  0.       0.      -0.70592  0.     ;
  0.91964 -0.5326   0.       0.     ;
 -0.35243  0.       0.      -1.70911;
  0.       0.       0.       0.62856;
 -1.89664  0.       0.       0.     ;
  0.       0.60997  1.61092  0.     ];
X = [[-0.31087  0.21058 -0.04888  0.1909 ]
 [-0.22324  0.09967  0.38207 -0.21312]
 [-0.13922 -0.23856  0.7192   0.17835]
 [ 0.09209 -0.049    0.13422 -0.44637]];
LearnHeat(X,20,0.1,0.1,[1;2],0)
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