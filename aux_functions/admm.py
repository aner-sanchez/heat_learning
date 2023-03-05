import numpy as np
import cvxpy as cp

def admm(X, Lt, gradient, Htp1, taut, dt, beta, verbose=False):
    S = len(taut)
    # variable
    L = cp.Variable(Lt.shape)
    N = Lt.shape[0]
    # constraints
    constraints = [cp.trace(L)== N, L.T==L, L@np.ones((N,1))==np.zeros((N,1))]
    for i in range(N-1):
        constraints += [L[i][i+1:]<=0]
    # objective
    obj = cp.Minimize(cp.trace(gradient.T@(L-Lt)) \
                      + dt/2*(cp.norm(L-Lt, 'fro')**2) + beta*(cp.norm(L, 'fro')**2))
    # solve problem
    prob = cp.Problem(obj, constraints)
    prob.solve(verbose=verbose, solver=cp.SCS, scale = 1000, use_indirect = False)
    if L.value is None:
        prob.solve(verbose=verbose, solver=cp.MOSEK)
    return L.value

L = admm(X, Lt, gradient, Htp1, taut, dt, beta)