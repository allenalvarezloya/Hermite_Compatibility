x = [1/5; 1/10; 1/20; 1/40; 1/80];
semilogy(x,SM2,x,SB2,x,SM3,x,SB3,'LineWidth',2)
set(gca,'FontSize',16)
xlabel('h','FontSize',16);
ylabel('Condition Number','FontSize',16);
legend('Non-eq Matrix m = 2','eq Matrix m = 2','Non-eq Matrix m = 3','eq Matrix m = 3','FontSize',12)
saveas(gcf,'polarConditioning','epsc')