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
H = zeros(N,1);
L1 = 0.0;
L2 = 0.0;
Linf = 0.0;
T = 0.0;

k31 = 6.3801618959239383506237;
k33 = 13.01520072169843441983;
gamma = k31/k33;
alpha = 1-gamma;

%% m = 2
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
            J = 2*pi*alpha*(alpha*r(i)+gamma);
            Uactual(i,j) = besselj(3,k33*(alpha*r(i)+gamma))*cos(6*pi*s(j))*cos(k33*T);
            L1 = L1 + J*abs(U(i,j) - Uactual(i,j));
            L2 = L2 + J*(U(i,j) - Uactual(i,j)).^2;
            if abs(U(i,j) - Uactual(i,j)) > Linf
                Linf = abs(U(i,j) - Uactual(i,j));
            end
        end
    end
    L1errors(k) = L1*hr*hs;
    L2errors(k) = sqrt(L2*hr*hs);
    Linferrors(k) = Linf;
    H(k) = 10*hr;
    L1 = 0.0;
    L2 = 0.0;
    Linf = 0.0;
end
L2_m2 = L2errors;

%% m = 3
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
            J = 2*pi*alpha*(alpha*r(i)+gamma);
            Uactual(i,j) = besselj(3,k33*(alpha*r(i)+gamma))*cos(6*pi*s(j))*cos(k33*T);
            L1 = L1 + J*abs(U(i,j) - Uactual(i,j));
            L2 = L2 + J*(U(i,j) - Uactual(i,j)).^2;
            if abs(U(i,j) - Uactual(i,j)) > Linf
                Linf = abs(U(i,j) - Uactual(i,j));
            end
        end
    end
    L1errors(k) = L1*hr*hs;
    L2errors(k) = sqrt(L2*hr*hs);
    Linferrors(k) = Linf;
    H(k) = 10*hr;
    L1 = 0.0;
    L2 = 0.0;
    Linf = 0.0;
end
L2_m3 = L2errors;
%% Plot
loglog(H,L2_m2,H,L2_m3,H,0.5*H.^4,'--k',H,20^(-2)*H.^6,'--k','LineWidth',2)
set(gca,'FontSize',16)
xlabel('h','FontSize',16);
ylabel('L2 Error','FontSize',16);
legend('m = 2','m = 3','','','FontSize',12)
saveas(gcf,'polarErrorConservative','epsc')