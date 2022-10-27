n = 3;
numGs = 1;


F = [-6, 4, -1; 4, -6, 5; -1, 5, -7]
G = [1, 0 , -2; 1, -1, 0; 0, 1, -1]

eig(F)
eig(G)



p = 0.1;
eps = 0.01;

seed = 2;
M = BMISolution(F, G, p, eps, seed); 