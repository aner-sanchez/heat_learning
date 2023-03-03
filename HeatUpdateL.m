% Implementation of proximal operator prox_ct^f(L(t)-grad(Z(L,H,tau),L)/dt)
% which is minimizer of <L-L(t),grad(Z(L,H,tau),L)>+dt/2
% *norm(L-L(t),fro)^2 + beta*norm(L,fro)^2 found in https://arxiv.org/abs/1611.01456
% page 5, update of L, we will use vech and vec to rewrite the problem in
% style found in

function L = HeatUpdateL(L_old,grad_L,lip_d,params)
    if exist('M','var')==0
        Mdup = duplication_matrix(params.n);
    end
    P = (lip_d/2+params.beta).*Mdup'*Mdup;
    % I don't know if this is good
    f = (grad_L(:))'*Mdup - lip_d*L_old(:)'*Mdup;
    if exist('A','var')||exist('B','var')==0
        A = eq_constraints_quad(params.n);
        B = ineq_constraints_quad(params.n);
    end
    options = optimset('Display','off');
    Ltemp = quadprog(P,f,B,zeros(size(B,1),1),A,[zeros(size(A,1)-1,1);params.n],...
        [],[],[],options);
    L = reshape(Mdup*Ltemp,params.n,params.n);
end

%%
function A = eq_constraints_quad(n)
    A = zeros(n+1,n*(n+1)/2);
    for i=1:n
        temp = zeros(n,n);
        temp(i,:) = ones(1,n); temp(:,i) = ones(n,1);
        A(i,:) = vech(temp);
    end
    % this is to lock trace to measure n, I have good results without this
    % and just normalizing the laplacian afterwards
    for i=1:n
        temp = zeros(n,n);
        temp(i,i) = 1;
        A(n+1,:) = A(n+1,:) + vech(temp)';
    end
end
%%
function B = ineq_constraints_quad(n)
    B = zeros(n*(n+1)/2-n,n*(n+1)/2);
    k = 2;
    iter = 1;
    flag = n-1;
    while flag > 0
        for i=1:flag
            B(iter,k) = 1;
            iter = iter + 1;
            k = k+1;
        end
        k = k+1;
        flag = flag - 1;
    end
end
%%
function D = duplication_matrix(n)
    % code from https://es.mathworks.com/matlabcentral/answers/473737-efficient-algorithm-for-a-duplication-matrix
    m   = n * (n + 1) / 2;
    nsq = n^2;
    D   = spalloc(nsq, m, nsq);
    row = 1;
    a   = 1;
    for i = 1:n
       b = i;
       for j = 0:i-2
          D(row + j, b) = 1;
          b = b + n - j - 1;
       end
       row = row + i - 1;
       
       for j = 0:n-i
          D(row + j, a + j) = 1;
       end
       row = row + n - i + 1;
       a   = a + n - i + 1;
    end
end
%%
function vech = vech(matA)
    % Function return elements from and below the main diagonal, then
    % stacking by column to have a vector of size K*(K+1)/2.
    % For example:  matA = [b11 b12 b13;b21 b22 b23;b31 b32 b33] is (3x3)
    %               vech(matA) returns: [b11 b12 b13 b22 b23 b33]' is (6x1)
    % Author: Binh T. Pham. Date: 2018-02-15.
    [M,N] = size(matA);
    if (M == N)
        vech  = [];
        for ii=1:M
            vech = [vech; matA(ii:end,ii)];
        end
    else
         error('Input must be a symmetric matrix.')
    end
end