% interbank plots - add corridor theta
line([ln_THETA(1) ln_THETA(end)], [delta_r delta_r],'Color','k','LineStyle',corridorstyle,'LineWidth',corridorwidth);
line([ln_THETA(1) ln_THETA(end)], [0 0],'Color','k','LineStyle',corridorstyle,'LineWidth',corridorwidth);
line([0 0], [0 delta_r],'Color','k','LineWidth',centerwidth,'LineStyle',centerstyle);
text(ln_THETA((N_theta-1)/2),delta_r-6,'$\mathbf{r^{w}-r^{m}}$','FontSize',16);
axis([theta0_bot theta0_top -5 delta_r+5]);
xlabel('$\mathbf{ln(\theta_0)}$');


ax = gca;
setNumYTicks(ax, delta_r, 5);
grid on;

function setNumYTicks(ax,delta_r, N)
  %  limits = ax.YLim;
    ax.YTick = linspace(0, delta_r, N);
end