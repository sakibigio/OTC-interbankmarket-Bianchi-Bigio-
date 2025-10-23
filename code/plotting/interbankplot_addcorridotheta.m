% interbank plots - add corridor theta
line([ln_THETA(1) ln_THETA(end)], [delta_r delta_r],'Color','k','LineStyle',corridorstyle,'LineWidth',corridorwidth);
line([ln_THETA(1) ln_THETA(end)], [0 0],'Color','k','LineStyle',corridorstyle,'LineWidth',corridorwidth);
line([0 0], [0 delta_r],'Color','k','LineWidth',centerwidth,'LineStyle',centerstyle);
text(ln_THETA((N_theta-1)/2),delta_r-6,'$\mathbf{r^{w}-r^{m}}$','FontSize',16);
