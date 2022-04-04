% This script computes the L1, L2 and Linf
% Errors and creates a table
u = dir('uplot*.ext');
x = dir('Yrefined*.ext');
y = dir('Xrefined*.ext');
N = length(u);
L1errors = zeros(N,1);
L2errors = zeros(N,1);
Linferrors = zeros(N,1);
L1 = 0.0;
L2 = 0.0;
Linf = 0.0;
hss = 0.1;
for k = 1:N
    U = load(u(k).name);
    X = load(x(1).name);
    Y = load(y(1).name);
    surf(X,Y,U)
end

