% implementing backtracking algorithm from https://github.com/semink/learnHeat
% mine doesn't work, just blows up / doesn't lower loss with L update

function L_new = backtrackingL(X,L,H,tau,params,lip,gamma)
    eta = 1.1; c2 = 0.01; k = 1; cond = false; iter = 1;
    grad_L = HeatGradL(X,L,H,tau,params);
    while cond==false && iter < 1000
        c2 = (eta^k)*c2;
        lip.d = gamma.d*c2; N = params.n;
%       Calling python is a pain, its really version dependant, I will try
%       to make it in matlab
%         Ltp1 = pyrunfile("admm.py","L",X=X, Lt=L, gradient=gradient, Htp1=H, ...
%             taut=tau, dt=lip.d, beta=params.beta);
        cvx_begin quiet;
            variable L_new(N, N);
            minimize(trace(grad_L'*(L_new - L)) +...
                (lip.d/2)*pow_pos(norm(L_new-L,"fro"),2) +...
                params.beta*pow_pos(norm(L_new,"fro"),2));
            subject to;
                trace(L_new) == N;
                L_new == L_new';
                for i=1:N
                    for j=1:N
                        if i~=j
                            L_new(i,j) <= 0;
                        end
                    end
                end
                L_new*ones(N,1) == zeros(N,1);
        cvx_end;
        next_L = L_new;
        k = k + 1;
        cond = L_condition(X,L,next_L,grad_L,H,tau,lip.d);
    end
    assert(iter<1000,"iterations reached 1000, something went wrong at backtrackingL.m");
end

%% AUXILIARY FUNCTIONS

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
