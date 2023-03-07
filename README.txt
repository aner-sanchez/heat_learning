This is code to emulate the algorithm exposed in https://arxiv.org/abs/1611.01456.
To start using it type HeatStart in the command line, if it doesn't work, then just initiate cvx as you can, it doesn't even need to be in the same folder, and type DemoHeat afterwards.

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

random_graph_RBF.m

Creates a random graph just as in the paper for testing purposes using a RBF.

randsignal.m

Creates a random sparse signal on a graph given the Laplacian, creates the dictionary
and linearly combines three random columns with standard gaussian coefficients
m times.

TesterHeat.m

Secondary script to test if everything is working. Can generate random inputs.

