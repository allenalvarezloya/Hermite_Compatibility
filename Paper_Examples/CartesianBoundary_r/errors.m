% This script computes the L1, L2 and Linf
% Errors and creates a table
u2 = dir('u2refined*.ext');
u3 = dir('u3refined*.ext');
x = dir('Xrefined*.ext');
y = dir('Yrefined*.ext');
N = length(u2);
L1errors = zeros(N,1);
L2errors = zeros(N,1);
Linferrors = zeros(N,1);
L1 = 0.0;
L2 = 0.0;
Linf = 0.0;
H = zeros(N,1);
for k = 1:N
    U = load(u2(k).name);
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
            Uactual(i,j) = sin(2*pi*X(i,j))*sin(2*pi*Y(i,j));
            L1 = L1 + abs(U(i,j) - sin(2*pi*X(i,j))*sin(2*pi*Y(i,j)));
            L2 = L2 + (U(i,j) - sin(2*pi*X(i,j))*sin(2*pi*Y(i,j))).^2;
            if abs(U(i,j) - sin(2*pi*X(i,j))*sin(2*pi*Y(i,j))) > Linf
                Linf = abs(U(i,j) - sin(2*pi*X(i,j))*sin(2*pi*Y(i,j)));
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
L2_m2 = L2errors;  
m = 2;
Ncells = 10*2.^(0:N-1)';

conv1 = log2(L1errors(1:end-1)./L1errors(2:end));
conv2 = log2(L2errors(1:end-1)./L2errors(2:end));
convinf = log2(Linferrors(1:end-1)./Linferrors(2:end));

%% Error Table
% fprintf("\\begin{table}[ht]\n ")
% fprintf("\\begin{center} \n ")
% fprintf("\\begin{tabular}{|c|c|c|c|c|c|c|} \n")
% fprintf("\\hline \n")
% fprintf(" n & $L_1$ error & Convergence & $L_2$ error & Convergence & $L_{\\infty}$ error & Convergence \\\\ \n ")
% fprintf("\\hline \n")
% fprintf(sprintf(' %i & %3.2e & - & %3.2e & - & %3.2e & - \\\\\\\\ \n',Ncells(1),L1errors(1),L2errors(1),Linferrors(1)))
% for z = 2:length(Ncells)
%     fprintf(sprintf('%i & %3.2e & %3.2f & %3.2e & %3.2f & %3.2e & %3.2f \\\\\\\\ \n',Ncells(z),L1errors(z),conv1(z-1),L2errors(z),conv2(z-1),Linferrors(z),convinf(z-1)))
% end
% fprintf("\\hline \n")
% fprintf("\\end{tabular} \n")
% fprintf("\\end{center} \n")
% fprintf("\\end{table} \n")

%% m = 3
L1errors = zeros(N,1);
L2errors = zeros(N,1);
Linferrors = zeros(N,1);
L1 = 0.0;
L2 = 0.0;
Linf = 0.0;
L2_actual = 0.0;
H = zeros(N,1);
for k = 1:N
    U = load(u3(k).name);
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
            Uactual(i,j) = sin(2*pi*X(i,j))*sin(2*pi*Y(i,j));
            L1 = L1 + abs(U(i,j) - sin(2*pi*X(i,j))*sin(2*pi*Y(i,j)));
            L2 = L2 + (U(i,j) - sin(2*pi*X(i,j))*sin(2*pi*Y(i,j))).^2;
            L2_actual = L2_actual + abs(Uactual(i,j));
            if abs(U(i,j) - sin(2*pi*X(i,j))*sin(2*pi*Y(i,j))) > Linf
                Linf = abs(U(i,j) - sin(2*pi*X(i,j))*sin(2*pi*Y(i,j)));
            end
        end
    end
    L1errors(k) = L1*hr*hs;
    L2errors(k) = sqrt(L2*hr*hs);
    L2errors_actual = sqrt(L2_actual*hr*hs);
    Linferrors(k) = Linf;
    H(k) = 10*hr;
    L1 = 0.0;
    L2 = 0.0;
    Linf = 0.0;
end


m = 3;
Ncells = 10*2.^(0:N-1)';

conv1 = log2(L1errors(1:end-1)./L1errors(2:end));
conv2 = log2(L2errors(1:end-1)./L2errors(2:end));
convinf = log2(Linferrors(1:end-1)./Linferrors(2:end));

%% Error Table
% fprintf("\\begin{table}[ht]\n ")
% fprintf("\\begin{center} \n ")
% fprintf("\\begin{tabular}{|c|c|c|c|c|c|c|} \n")
% fprintf("\\hline \n")
% fprintf(" n & $L_1$ error & Convergence & $L_2$ error & Convergence & $L_{\\infty}$ error & Convergence \\\\ \n ")
% fprintf("\\hline \n")
% fprintf(sprintf(' %i & %3.2e & - & %3.2e & - & %3.2e & - \\\\\\\\ \n',Ncells(1),L1errors(1),L2errors(1),Linferrors(1)))
% for z = 2:length(Ncells)
%     fprintf(sprintf('%i & %3.2e & %3.2f & %3.2e & %3.2f & %3.2e & %3.2f \\\\\\\\ \n',Ncells(z),L1errors(z),conv1(z-1),L2errors(z),conv2(z-1),Linferrors(z),convinf(z-1)))
% end
% fprintf("\\hline \n")
% fprintf("\\end{tabular} \n")
% fprintf("\\end{center} \n")
% fprintf("\\end{table} \n")

L2_m3 = L2errors;  

loglog(H,L2_m2,H,L2_m3,H,8*H.^3,'--k',H,2.5*H.^5,'--k','LineWidth',2)
set(gca,'FontSize',16)
xlabel('h','FontSize',16);
ylabel('L2 Error','FontSize',16);
legend('m = 2','m = 3','','','FontSize',12)
saveas(gcf,'cartesianErrorDissipative','epsc')