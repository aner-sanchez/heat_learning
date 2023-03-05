clc, clear all
disp("Let's start by creating a random graph using a RBF"); pause(2);
n = input("How many vertices do you want the graph to have?\n" + ...
    "Write and press enter: ");
L = randgraph_RBF(n);
W = -(L-diag(diag(L)));
disp("I will plot the graph at the end of the demostration")
disp("Or you can just type plot(G) in the end"); pause(1);
G = graph(W); plot(G);
disp("As you will see in the Workspace we have the following variables");
disp("G: the graph object")
disp("L: the Laplacian of the graph")
disp("n: the number of vertices")
disp("W: the adjacency matrix")
disp("Now I'm going to choose for you something " + ...
    "we are going to use the following tau vector")
tau = [2.5;4]
disp("that is: there are two hidden diffusion processes over the graph" + ...
    "and those are the diffusion constants")
disp("Let's create m random sparse signals over the graph")
m = input("Choose m: ");
X = randsignal(L,m,tau);
clc
disp("We've just created a random signal, look at it 5 seconds")
X
pause(5)
disp("suppose we only have the signal X, how can we know what was the Laplacian?")
disp("not only that, but the rate of diffusions and the sparse signal")
pause(2)
disp("LearnHeat.m will attempt to recover, input the following parameters")
alpha = input("alpha near zero (example: 0.1) (to control sparsity): ")
beta = input("beta near zero (example: 0.1) (to control sparsity too): ")
iter = input("iter: number of iterations (you can put 20) <- ")
quiet = input("type 0 if you want information through the optimization process (you surely do) ")
disp("I will choose a tau for you, which is not the real tau")
tau2 = [1;2]
disp("now is up to you: type the following in the terminal")
disp("LearnHeat(X,iter,alpha,beta,tau2,quiet)")