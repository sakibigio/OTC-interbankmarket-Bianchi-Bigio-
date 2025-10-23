line([round(1) round(end)], [delta_r delta_r],'Color','k','LineStyle',string(corridorstyle),'LineWidth',corridorwidth);
line([round(1) round(end)], [0 0],'Color','k','LineStyle',corridorstyle,'LineWidth',corridorwidth);
line([0 0], [0 delta_r],'Color','k','LineWidth',centerwidth,'LineStyle',centerstyle);
text(round((N_theta-1)/2),delta_r-6,'$\mathbf{r^{w}-r^{m}}$','FontSize',16);
axis([round(1) round(end) -5 delta_r+5]);

ax = gca;
setNumYTicks(ax, delta_r, 5);
grid on;

function setNumYTicks(ax,delta_r, N)
  %  limits = ax.YLim;
    ax.YTick = linspace(0, delta_r, N);
end