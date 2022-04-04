loglog(H,L2_m2,H,L2_m3,H,10e3*H.^4,'--k',H,20e4*H.^6,'--k','LineWidth',2)
set(gca,'FontSize',16)
xlabel('h','FontSize',16);
ylabel('L2 Error','FontSize',16);
legend('m = 2','m = 3','','','FontSize',12)
saveas(gcf,'cartesianErrorConservative','epsc')