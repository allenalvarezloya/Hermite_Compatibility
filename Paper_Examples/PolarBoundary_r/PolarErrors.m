loglog(H(1:8:end),L2_m2(1:8:end),H(1:8:end),L2_m3(1:8:end),H(1:8:end),8e4*H(1:8:end).^5,'--k',H(1:8:end),3.5e3*H(1:8:end).^3,'--k','LineWidth',2)
set(gca,'FontSize',16)
xlabel('h','FontSize',16);
ylabel('L2 Error','FontSize',16);
legend('m = 2','m = 3','','','FontSize',12)
saveas(gcf,'polarErrorDissipative','epsc')