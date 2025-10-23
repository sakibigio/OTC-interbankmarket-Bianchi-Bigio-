% interbankplot_addstaticrate
text(ln_THETA(80),(1-matchtech.eta)*delta_r-6,'$(1-\eta)(\mathbf{r^{w}-r^{m}})$','FontSize',16);
line([ln_THETA(1) ln_THETA(end)], (1-matchtech.eta)*[delta_r delta_r],'Color','k','LineStyle','-','LineWidth',1);
