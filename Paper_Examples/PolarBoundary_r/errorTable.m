% This script computes the L1, L2 and Linf
% Errors and creates a table
u = dir('urefined*.ext');
x = dir('Xrefined*.ext');
y = dir('Yrefined*.ext');
N = length(u);
L1errors = zeros(N,1);
L2errors = zeros(N,1);
Linferrors = zeros(N,1);
H = zeros(N,1);
L1 = 0.0;
L2 = 0.0;
Linf = 0.0;
hss = 0.1;

k31 = 6.3801618959239383506237;
k33 = 13.01520072169843441983;
gamma = k31/k33;
alpha = 1-gamma;
for k = 1:N
    U = load(u(k).name);
    X = load(x(k).name);
    Y = load(y(k).name);
    Uactual = 0*U;
    S = size(U);
    nr = S(1);
    ns = S(2);
    hr = 1.0/(nr-1);
    hs = 1.0/(ns-1);
    r = (0:nr-1)*hr;
    s = (0:ns-1)*hs;
    for j = 1:ns
        for i = 1:nr
            J = 2*pi*alpha*(alpha*r(i)+gamma);
            Uactual(i,j) = besselj(3,k33*(alpha*r(i)+gamma))*cos(6*pi*s(j));
            L1 = L1 + J*abs(U(i,j) - besselj(3,k33*(alpha*r(i)+gamma))*cos(6*pi*s(j)));
            L2 = L2 + J*(U(i,j) - besselj(3,k33*(alpha*r(i)+gamma))*cos(6*pi*s(j))).^2;
            if abs(U(i,j) - besselj(3,k33*(alpha*r(i)+gamma))*cos(6*pi*s(j))) > Linf
                Linf = abs(U(i,j) - besselj(3,k33*(alpha*r(i)+gamma))*cos(6*pi*s(j)));
            end
        end
    end
    L1errors(k) = L1*hr*hs;
    L2errors(k) = sqrt(L2*hr*hs);
    Linferrors(k) = Linf;
    H(k) = hr;
    L1 = 0.0;
    L2 = 0.0;
    Linf = 0.0;
end
plot(X(1:10:end,1:10:end),Y(1:10:end,1:10:end),'k',X(1:10:end,1:10:end)',Y(1:10:end,1:10:end)','k')

m = 2;
Ncells = 5*2.^(0:N-1)';

conv1 = log2(L1errors(1:end-1)./L1errors(2:end));
conv2 = log2(L2errors(1:end-1)./L2errors(2:end));
convinf = log2(Linferrors(1:end-1)./Linferrors(2:end));

%% Error Table
fprintf("\\begin{table}[ht]\n ")
fprintf("\\begin{center} \n ")
fprintf("\\begin{tabular}{|c|c|c|c|c|c|c|} \n")
fprintf("\\hline \n")
fprintf(" n & $L_1$ error & Convergence & $L_2$ error & Convergence & $L_{\\infty}$ error & Convergence \\\\ \n ")
fprintf("\\hline \n")
fprintf(sprintf(' %i & %3.2e & - & %3.2e & - & %3.2e & - \\\\\\\\ \n',Ncells(1),L1errors(1),L2errors(1),Linferrors(1)))
for z = 2:length(Ncells)
    fprintf(sprintf('%i & %3.2e & %3.2f & %3.2e & %3.2f & %3.2e & %3.2f \\\\\\\\ \n',Ncells(z),L1errors(z),conv1(z-1),L2errors(z),conv2(z-1),Linferrors(z),convinf(z-1)))
end
fprintf("\\hline \n")
fprintf("\\end{tabular} \n")
fprintf("\\end{center} \n")
fprintf("\\end{table} \n")