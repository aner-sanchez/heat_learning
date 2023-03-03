This is code to emulate the algorithm exposed in https://arxiv.org/abs/1611.01456.

HeatDict.m

Computes the heat dictionary [exp(-tau_i(L))]_i=1:s
the exponential of minus tau laplacian for each underlying diffusion process.

HeatGradL.m

Computes the gradient of the objective function w.r.t L.

HeatGradtau.m

Computes the gradient of the objective function w.r.t tau. There is no HeatGradH.m
because it is automatic.

HeatLipL.m

Computes Lipschitz constant for the L-step using a backtracking algorithm.

HeatLiptau.m

Computes Lipschitz constant for the tau-step.

HeatUpdateL.m

Computes the update of L given the gradient and the lipschitz constant, it is some
sort of proximal operator, it is a quadratic program using vectorization of matrices.

LearnHeat.m

Main script, collects everything and alternates between minimization of L, H,
tau.

TesterHeat.m

Secondary script to test if everything is working. Can generate random inputs.

