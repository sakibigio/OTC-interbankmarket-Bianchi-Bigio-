% interbank plots - add corridor theta
line([LAMBDA_vec(1) LAMBDA_vec(end)], [delta_r delta_r],'Color','k','LineWidth',corridorwidth,'LineStyle',corridorstyle);
line([LAMBDA_vec(1) LAMBDA_vec(end)], [0 0],'Color','k','LineWidth',corridorwidth,'LineStyle',corridorstyle);
text(LAMBDA_vec(5),delta_r-5,'$\mathbf{r^{w}-r^{m}}$','FontSize',16);
axis([LAMBDA_vec(1) LAMBDA_vec(end) -5 delta_r+5]);
xlabel('$\mathbf{\lambda}$')

ax = gca;
setNumYTicks(ax, delta_r, 5);
grid on;

function setNumYTicks(ax,delta_r, N)
  %  limits = ax.YLim;
    ax.YTick = linspace(0, delta_r, N);
end