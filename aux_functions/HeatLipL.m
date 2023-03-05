% https://arxiv.org/abs/1611.01456
% Computing Lipschitz constant for L, C_2(H,tau) by backtracking, see
% appendix B-B
% Input: eta=1.1 (prevalue), initial guess for C_2(H,tau), k=1
% Output: estimate Lipschitz constant C_2(H,tau)

function [L_new,lip_d] = HeatLipL(X,L_old,L_new,grad_L,H,tau,params,gamma,lip_d,eta)
    if nargin < 10; eta = 1.1; end
    if nargin < 9; lip_d = 1; end
    k = 1; iter = 1; C_2 = lip_d;
    %% while (23) is false, do
%     assert(all(size(lip_d)==[1,1]),['The lipschitz constant being passed to HeatLipL as an initial guess...' ...
%         'is not a scalar'])
    while L_condition(X,L_old,L_new,grad_L,H,tau,C_2) == 0 && iter < 1e3
        C_2 = (eta.^k)*C_2; k = k + 1; iter = iter + 1; lip_d = gamma.d*C_2;
        L_new = HeatUpdateL(L_old,grad_L,lip_d,params);
    end
    if iter == 1e3
        warning("Max iterations were reached in HeatLipL.m");
    end
    % we multiply by gamma > 1 here as iterations are finished and this
    % must be to ensure being strictly above the lipschitz constant
    lip_d = gamma.d*lip_d;
    %assert(all(size(lip_d)==[1,1]),'Exiting lip_d from HeatLipL.m is not a scalar')
end

function L_bool = L_condition(X,L_old,L_new,grad_L,H,tau,lip_d)
    % this is tricky, of course everything should be a scalar, this is my
    % interpretation of A'A as vec(A)'vec(A)
    % my option
%     L_bool = (Z(X,L_new,H,tau) <= Z(X,L_old,H,tau) + (grad_L(:))'*(L_new(:)-L_old(:)) ...
%         + (lip_d/2)*norm(L_new-L_old,"fro").^2);
    % EPFL option
    lhs = Z(X,L_new,H,tau);
    rhs = matrix_inner(grad_L,L_new-L_old)+(lip_d/2)*norm(L_new-L_old,"fro").^2;
    L_bool = (lhs <= rhs); %§
end

%% Distance from current signal estimation and dictionary to real noisy signal
% Page 4

function Z = Z(X,L,H,tau)
    D = HeatDict(L,tau);
    Z = norm(X - D*H,"fro").^2;
end

function inner = matrix_inner(A,B)
    % I think there is another way of writing this as vec(B)^T * vec(A)
    inner = trace(B'*A);
end