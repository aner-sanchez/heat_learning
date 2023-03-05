% Learns a graph from signals assuming an underlying diffusion: see
% Algorithm 1 in https://arxiv.org/abs/1611.01456
% Input: Signal set X
% Optional inputs:
    % number of iterations iter
    % parameters alpha and beta
    % vector of diffusion processes tau
    % Lipschitz coefficients augmentation gamma
% Output: Sparse signal representations H, graph Laplacian L, diffusion
% parameter tau
% IMPORTANT: algorithm really depends on Lipschitz constants, if they are
% not set up correctly output wont make any sense

function [L,H,tau] = LearnHeat(X,iter,alpha,beta,tau,quiet,gamma)
    %% params
    aux = 1.1; if nargin < 7; gamma = struct('c',aux,'d',aux,'e',aux); end
    if nargin < 6; quiet = 1; end
    if nargin < 5; tau = [2.5;4]; end
    if nargin < 4; beta = 0.01; end
    if nargin < 3; alpha = 0.1; end
    if nargin < 2; iter = 20; end
    if nargin == 0
        error('Signal slot empty');
    end
    params = struct('iter',iter,'alpha',alpha,'beta',beta,'tau',tau,...
        'n',size(X,1),'S',size(tau,1),'m',size(X,2));
    %% Initialization
    L = zeros(params.n); D = HeatDict(L,tau); H = zeros(params.n.*params.S,params.m);
    % store lipschitz constants in a struct
    lip = struct('c',0,'d',0,'e',0); lip.d = 0.1; flag = 1;
    %% Iterations, alternating between H, L, tau
    for t=1:iter
        fprintf("Iteration Nº%d\n",t);
        %% save L,H,tau just in case cost goes up
        if flag ~=0
            savL = L; savH = H; savtau = tau;
            savcost = HeatCalCost(X,L,H,tau,params);
        end
        %% H
        lip.c = gamma.c.*norm(2.*(D')*D,"fro"); %§
        % solve problem (8) to update H
        grad_H = -2*(D')*(X-D*H);
        prox_ct_f = @(z) sign(z).*(max(abs(z)-alpha/lip.c,0));
        H = prox_ct_f(H - grad_H./lip.c);
        if quiet == 0; fprintf("H updated, cost: %f\n",HeatCalCost(X,L,H,tau,params)); end
        %% L, most complex part of the algorithm
        % we should put here a better L updater, 
        L = backtrackingL(X,L,H,tau,params,lip,gamma);
        if quiet == 0; fprintf("L updated, cost: %f\n",HeatCalCost(X,L,H,tau,params)); end
        %% tau
        lip.e = gamma.e.*HeatLiptau(X,L,H,params); % this makes the algorithm not work if you dont put anything correct
        grad_tau = HeatGradtau(X,L,H,tau,params);
        tau = max(-(grad_tau-lip.e*tau)/lip.e,0);
        % update dictionary
        D = HeatDict(L,tau,D);
        if quiet == 0; fprintf("Tau updated, cost: %f\n",HeatCalCost(X,L,H,tau,params)); end
        %% print current cost
        current_cost = HeatCalCost(X,L,H,tau,params);
        if current_cost > savcost
            fprintf("Warning: cost goes up with iterations");
            flag = 0;
        else
            flag = 1;
        end
        fprintf("\n\n");
    end
    if flag == 0
        L = savL; H = savH; tau = savtau;
        warning("We were not able to descent monotonically.");
    end
    L(L<1e-4) = 0; H(H<1e-4) = 0;
end