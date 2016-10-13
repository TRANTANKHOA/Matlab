%% Plotting whole trajectory of alpha
alpha_pct = quantile( alpha_M, [0.95, 0.05]);
f = figure('name',window);
plot(alpha_pct')
hold on;
plot(alpha,'-.')
legend('95 %','5 %',['x_{' num2str(t) '}' ],'Location','southeast')
xlabel('Time');ylabel(['$x_{' num2str(t) '}$' ])
title(['Posterior density of $x_{' num2str(t) '}$'])

