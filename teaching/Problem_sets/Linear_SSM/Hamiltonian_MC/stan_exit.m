    %% Printing notification
    fprintf('\n');
    beep;
    if any(fit.exit_value)
        fprintf('Stan failed with exitValue = \n');
        disp(fit.exit_value)
    else
        fprintf('Stan succeeded. \n');
    end
    fprintf('\n');

    %% Plotting sigma
    window = 'Stan output';
    plot_time_series(fit.extract('permuted',true).logsigma,'$\log(\sigma^2)$',10,window);
    title('Sampling $\theta$ using STAN')
    %% Plotting tau
    window = 'Stan output';
    plot_time_series(fit.extract('permuted',true).logtau,'$\log(\tau^2)$',10,window);
    title('Sampling $\theta$ using STAN')
    %% Plotting phi
    plot_time_series(fit.extract('permuted',true).phi,'$\phi$',10,window);
    title('Sampling $\phi$ using STAN')
    %% Plotting generated data y
    Y = fit.extract('permuted',true).y_tilde;
    plot_y_T
    %% Plotting trace of x_t
t = randi(sys.T);
x_M = fit.extract('permuted',true).x;
plot_time_series(x_M(:,t),['$x_{' num2str(t) '}$' ],10,window);
%% Plotting whole trajectory of x_T
x_pct = quantile( x_M, [0.95, 0.05]);
figure
plot(x_pct')
hold on;
plot(sys.x,'-.')
legend('95 %','5 %','x_t')
xlabel('Time');ylabel('$x$ - Series')
title('State trajectory $x_t$')
