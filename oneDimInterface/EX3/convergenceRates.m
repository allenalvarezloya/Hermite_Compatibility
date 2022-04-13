close all

U1 = dir('u1m1oversampled*.ext');
U2 = dir('u2m1oversampled*.ext');
N1 = length(U1);
L2_m1 = zeros(N1,1);
C = 1.0;
D = 2.0;
T = 2*C/(D+C);
R = T-1;
L2_m2 = zeros(N1,1);
H = zeros(N1,1);
for k = 1:N1
    u1 = load(U1(k).name);
    n1 = length(u1)-1;
    h1 = 2.0/n1;
    x1 = -1.0 + h1*(0:n1)';
    u2 = load(U2(k).name);
    n2 = length(u2)-1;
    h2 = 2.0/n2;
    x2 = -1.0 + h2*(0:n2)';
    u1 = u1.*(x1<=0.0);
    u2 = u2.*(x2>=0.0);
    u2ex = T*0.5*exp(-500*((x2-D*0.3)/D).^2).*(x2>=0.0);
    u1ex = R*0.5*exp(-500*((x1+C*0.3)/C).^2).*(x1<=0.0) - 0.5*exp(-500*((x1+C*0.7)/C).^2).*(x1<=0.0);
    L2_m2(k) = sqrt(h1*sum((u1 - u1ex).^2)) + sqrt(h2*sum((u2 - u2ex).^2));
    H(k) = h1;
end



U1 = dir('u1m2oversampled*.ext');
U2 = dir('u2m2oversampled*.ext');
N1 = length(U1);
L2_m2 = zeros(N1,1);
C = 1.0;
D = 2.0;
T = 2*C/(D+C);
R = T-1;
L2_m2 = zeros(N1,1);
H = zeros(N1,1);
for k = 1:N1
    u1 = load(U1(k).name);
    n1 = length(u1)-1;
    h1 = 2.0/n1;
    x1 = -1.0 + h1*(0:n1)';
    u2 = load(U2(k).name);
    n2 = length(u2)-1;
    h2 = 2.0/n2;
    x2 = -1.0 + h2*(0:n2)';
    u1 = u1.*(x1<=0.0);
    u2 = u2.*(x2>=0.0);
    u2ex = T*0.5*exp(-500*((x2-D*0.3)/D).^2).*(x2>=0.0);
    u1ex = R*0.5*exp(-500*((x1+C*0.3)/C).^2).*(x1<=0.0) - 0.5*exp(-500*((x1+C*0.7)/C).^2).*(x1<=0.0);
%     u2ex = 0.5*exp(-500*((x2-C*0.5)/C).^2).*(x2>=0.0);
%     u1ex = -0.5*exp(-500*((x1+C*0.5)/C).^2).*(x1<=0.0);
    L2_m2(k) = sqrt(h1*sum((u1 - u1ex).^2)) + sqrt(h2*sum((u2 - u2ex).^2));
    H(k) = h1;
end

U1 = dir('u1m3oversampled*.ext');
U2 = dir('u2m3oversampled*.ext');
N1 = length(U1);
L2_1 = zeros(N1,1);
C = 1.0;
D = 2.0;
T = 2*C/(D+C);
R = T-1;
L2_m3 = zeros(N1,1);
H = zeros(N1,1);
for k = 1:N1
    u1 = load(U1(k).name);
    n1 = length(u1)-1;
    h1 = 2.0/n1;
    x1 = -1.0 + h1*(0:n1)';
    u2 = load(U2(k).name);
    n2 = length(u2)-1;
    h2 = 2.0/n2;
    x2 = -1.0 + h2*(0:n2)';
    u1 = u1.*(x1<=0.0);
    u2 = u2.*(x2>=0.0);
    u2ex = T*0.5*exp(-500*((x2-D*0.3)/D).^2).*(x2>=0.0);
    u1ex = R*0.5*exp(-500*((x1+C*0.3)/C).^2).*(x1<=0.0) - 0.5*exp(-500*((x1+C*0.7)/C).^2).*(x1<=0.0);
%     u2ex = 0.5*exp(-500*((x2-C*0.5)/C).^2).*(x2>=0.0);
%     u1ex = -0.5*exp(-500*((x1+C*0.5)/C).^2).*(x1<=0.0);
    L2_m3(k) = sqrt(h1*sum((u1 - u1ex).^2)) + sqrt(h2*sum((u2 - u2ex).^2));
    H(k) = h1;
end

loglog(H,L2_m2,H,L2_m3,'LineWidth',2)
hold on
loglog(H,5e5*H.^3,'k--','LineWidth',2)
loglog(H,3e9*H.^5,'k--','LineWidth',2)
hold off
set(gca,'FontSize',20)
legend('m = 2 Dissipative','m = 3 Dissipative','','')
saveas(gcf,'travelingWaveDissipative','epsc')
