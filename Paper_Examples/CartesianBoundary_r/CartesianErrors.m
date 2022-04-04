loglog(H(1:3:end),L2_m2(1:3:end),H(1:3:end),L2_m3(1:3:end),H(1:3:end),25e4*H(1:3:end).^5,'--k',H(1:3:end),9e3*H(1:3:end).^3,'--k','LineWidth',2)
set(gca,'FontSize',16)
xlabel('h','FontSize',16);
ylabel('L2 Error','FontSize',16);
legend('m = 2','m = 3','','','FontSize',12)
saveas(gcf,'cartesianErrorDissipative','epsc')