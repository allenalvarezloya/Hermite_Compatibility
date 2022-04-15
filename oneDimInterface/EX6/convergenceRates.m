close all
U1 = dir('um1oversampled*.ext');
U2 = dir('um2oversampled*.ext');
U3 = dir('um3oversampled*.ext');
N = length(U1);
L2_m1 = zeros(N,1);
L2_m2 = zeros(N,1);
L2_m3 = zeros(N,1);
H = zeros(N,1);
for k = 1:N
    u1 = load(U1(k).name);
    u2 = load(U2(k).name);
    u3 = load(U3(k).name);
    n = length(u1)-1;
    h = 1.0/n;
    x = 0.0 + h*(0:n)';
    uex = sin(2*pi*(x-1.3));
    L2_m1(k) = sqrt(h*sum((u1(1:end-1) - uex(1:end-1)).^2));
    L2_m2(k) = sqrt(h*sum((u2(1:end-1) - uex(1:end-1)).^2));
    L2_m3(k) = sqrt(h*sum((u3(1:end-1) - uex(1:end-1)).^2));
    H(k) = h;
end
loglog(H,L2_m1,H,L2_m2,H,L2_m3,'LineWidth',2)
hold on 
loglog(H,22*H,'k--',H,22e2*H.^3,'k--',H,6e4*H.^5,'k--','LineWidth',2)
hold off
legend('m = 1 Dissipative','m = 2 Dissipative','m = 3 Dissipative','','','','FontSize',16)
set(gca,'FontSize',20)
xlabel('h','FontSize',20)
ylabel('L2 error','FontSize',20)
saveas(gcf,'dissipativeTraveling','epsc')