% compute gradient of Z(L,H,tau) w.r.t L as explained in https://arxiv.org/abs/1611.01456
% Appendix A-B
% Inputs: signal X, current L, H, tau
% Optional inputs:
    % grad_L to not reallocate again

%% APPENDIX A-B
function grad_L = HeatGradL(X,L,H,tau,params)
    % calloc grad_L
    grad_L = zeros(params.n);
    for s=1:params.S
        % size(H) = (n*S,m), size(X') = (m,n) so we need to take (n,m) from
        % H each step of s
        A = H(1+(s-1)*params.n:s*params.n,:)*X';
        % disp("size A"); size(A)
        % eigenvalue decomposition of exponent of expm
        grad_L = grad_L - 2.*grad_trace(A,-tau(s).*L,params);
        for s2=1:params.S
            % disp("we reached here")
            A = H(1+(s2-1)*params.n:s2*params.n,:)*(H(1+(s-1)*params.n:s*params.n,:)');
            grad_L = grad_L + grad_trace(A,-(tau(s)+tau(s2)).*L,params);
        end
    end
end

%% AUXILIARY FUNCTIONS
% following function computes gradient \nabla_L tr(A*expm(L));
function grad_tr = grad_trace(A,L,params)
    [chi,lambda] = eig(L);
    % Appendix A-B definition (16) of matrix B
    B = expm(lambda);
    for i=params.n
        for j=params.n
            if i~=j
                B(i,j) = (exp(lambda(i,i)-exp(lambda(j,j))))...
                    ./(lambda(i,i)-lambda(j,j));
            end
        end
    end
    % disp("size chi"); size(chi)
    assert(size(chi,1)==size(A,1));
    grad_tr = chi*((chi'*A*chi).*B)*chi';
end