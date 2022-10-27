n = 2;
numGs = 2;

a1 = 1;
a2 = 2;
b = 2;
s1 = 1;
s2 = 2;

F = [a1, -b; b, a2];
G = zeros(n, n, numGs);
G(:,:,1) = [s1, 0; 0, s1];
G(:,:,2) = [s2, 0; 0, s2];
eig(F)


p = 0.2;
eps = 0.001;

seed = 5;
M = BMISolution(F, G, p, eps, seed); 