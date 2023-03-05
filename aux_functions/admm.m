% trying to imitate admm.py with matlab

function L = admm(L_old, grad_L, lip_d, params)
    cvx_begin quiet
        variable L_new(params.n, params.n);
        minimize(trace(grad_L'*(L_new - L_old)) +...
            (lip_d/2)*pow_pos(norm(L_new-L_old,"fro"),2) +...
            params.beta*pow_pos(norm(L_new,"fro"),2));
        subject to;
            trace(L_new) == params.n;
            L_new == L_new';
            for i=1:params.n
                for j=1.params.n
                    L_new(i,j) <= 0;
                end
            end
            L_new*ones(params.n,1) == zeros(params.n,1);
    cvx_end
end