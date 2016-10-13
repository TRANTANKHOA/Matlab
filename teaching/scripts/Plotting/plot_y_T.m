%% Ploting y
Y_pct = quantile( Y, [0.95, 0.05]);
f = figure('name',window);
plot(Y_pct')
hold on;
plot(sys.y,'-.')
legend('95 %','5 %','y_t','Location','northeast')
xlabel('Time');ylabel('$y_t$')
title('Posterior predictive density of $y_t$')
axis auto
