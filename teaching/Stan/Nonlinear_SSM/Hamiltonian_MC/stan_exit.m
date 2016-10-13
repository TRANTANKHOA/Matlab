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

    %% Plotting theta
    window = 'HMC';
    plot_time_series(fit.extract('permuted',true).theta,'$\theta$',10,window);
    %% Plotting theta
    plot_time_series(fit.extract('permuted',true).phi,'$\phi$',10,window);
    title('Sampling $\phi$ using STAN')
    %% Plotting generated data y
    Y = fit.extract('permuted',true).y_tilde;
    plot_y_T
    title('Sampling $\theta$ using STAN')
    %% Plotting trace of alpha_t
t = randi(T);
alpha_M = fit.extract('permuted',true).alpha;
plot_time_series(alpha_M(:,t),['$\alpha_{' num2str(t) '}$' ],10,window);
title('Sampling $(\theta,\alpha_t)$ using STAN')
%% Plotting whole trajectory of alpha_T
alpha_pct = quantile( alpha_M, [0.95, 0.05]);
figure
plot(alpha_pct')
hold on;
plot(alpha,'-.')
legend('95 %','5 %','\alpha_t')
xlabel('Time');ylabel('$\alpha$ - Series')
title('Sampling $(\theta,\alpha_t)$ using STAN')
%% Cleaning up working folder
    delete anon_model.stan anon_model.hpp
