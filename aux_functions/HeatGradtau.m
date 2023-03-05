%% APPENDIX A-C
function grad_tau = HeatGradtau(X,L,H,tau,params)
    grad_tau = zeros(size(tau));
    for s=1:params.S
        grad_tau(s) = grad_tau(s) + 2*trace(H(1+(s-1)*params.n:s*params.n,:)*X'...
            *L*expm(-tau(s)*L));
        for s2=1:params.S
            grad_tau(s) = grad_tau(s) - 2*trace(H(1+(s2-1)*params.n:s2*params.n,:)...
                *(H(1+(s-1)*params.n:s*params.n,:)')*L*expm(-(tau(s)+tau(s2))*L));
        end
    end
end