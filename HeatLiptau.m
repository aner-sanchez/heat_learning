% https://arxiv.org/abs/1611.01456
% Computing Lipschitz constant for tau, appendix B-C, not trivial at all

function lip_e = HeatLiptau(X,L,H,params)
    % it is given by a max, this is O(S^2) complex
    running_max = 0; X_norm = norm(X,"fro");
    for s=1:params.S
        H_norm = norm(H(1+(s-1)*params.n:s*params.n,:),"fro");
        current_max = 2*H_norm*X_norm;
        for s2=1:params.S
            current_max = current_max + 4*H_norm*...
                norm(H(1+(s2-1)*params.n:s2*params.n,:),"fro");
        end
        running_max = max(running_max,current_max);
    end
    lip_e = norm(L).^2*running_max;
end