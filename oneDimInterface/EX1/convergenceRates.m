close all
U1 = dir('u1m1oversampled*.ext');
U2 = dir('u1m2oversampled*.ext');
U3 = dir('u1m3oversampled*.ext');
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
    uex = sin(2*pi*x)*cos(2*pi*1.3);
    L2_m1(k) = sqrt(h*sum((u1 - uex).^2));
    L2_m2(k) = sqrt(h*sum((u2 - uex).^2));
    L2_m3(k) = sqrt(h*sum((u3 - uex).^2));
    H(k) = h;
end
loglog(H,L2_m1,H,L2_m2,H,L2_m3,'LineWidth',2)
hold on 
loglog(H,50e1*H.^2,'k--',H,15e3*H.^4,'k--',H,30e4*H.^6,'k--','LineWidth',2)
hold off
legend('m = 1 Conservative','m = 2 Conservative','m = 3 Conservative','','','','FontSize',16)
set(gca,'FontSize',20)
xlabel('h','FontSize',20)
ylabel('L2 error','FontSize',20)
saveas(gcf,'conservative','epsc')
