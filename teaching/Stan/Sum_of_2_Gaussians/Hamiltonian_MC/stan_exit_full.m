%% Printing notification
fprintf('\n');
beep;
if any(fit_full.exit_value)
    fprintf('Stan failed with exitValue = \n');
    disp(fit_full.exit_value)
else
    fprintf('Stan succeeded. \n');
end
fprintf('\n');
%% Plotting theta
window = 'Hamiltonian Monte Carlo';
f = plot_time_series(fit_full.extract('permuted',true).theta,'$\upsilon$',10,window);
subplot(2,2,1)
set(gca,'Ylim',YL_theta)
subplot(2,2,2)
set(gca,'Ylim',YL_theta)
legend('off')
%% Printing
    w = 0.5; h = 0.35; w=round(w*1050); h=round(h*800); r=300; sz = [-1000 1000 w h];
    drawnow; undock(f); set(f,'Units','pixels','Position',sz);
    filename = '/Users/Eric/Dropbox/PhD/Reports/Literature Review/Pics/sum_2_Gauss_theta_HMC_full.pdf';
    export_fig(filename,f,'-a1',['-r',num2str(r)]);%, '-grey','-nocrop'
    saveas(f,[filename(1:end-4) '.fig'])

%% Plotting trace of alpha_t
alpha_M = fit_full.extract('permuted',true).alpha;
f = plot_time_series(alpha_M(:,t),['$x_{' num2str(t) '}$' ],10,window);
subplot(2,2,1)
set(gca,'Ylim',YL_alpha)
subplot(2,2,2)
set(gca,'Ylim',YL_alpha)
legend('off')

    %% Printing
    w = 0.5; h = 0.35; w=round(w*1050); h=round(h*800); r=300; sz = [-1000 1000 w h];
    drawnow; undock(f); set(f,'Units','pixels','Position',sz);
    filename = '/Users/Eric/Dropbox/PhD/Reports/Literature Review/Pics/sum_2_Gauss_alpha_HMC.pdf';
    export_fig(filename,f,'-a1',['-r',num2str(r)]);%, '-grey','-nocrop'
    saveas(f,[filename(1:end-4) '.fig'])
%% Plotting whole trajectory of alpha_T
plot_alpha_T
set(gca,'Ylim',YL_alpha_M)
%% Plotting generated data y
Y = fit_full.extract('permuted',true).y_tilde;
plot_y_T
legend('Location','northeastoutside')
set(gca,'Ylim',YL_y)
    %% Printing
    w = 1; h = 0.25; w=round(w*1050); h=round(h*800); r=300; sz = [-1000 1000 w h];
    drawnow; undock(f); set(f,'Units','pixels','Position',sz);
    filename = '/Users/Eric/Dropbox/PhD/Reports/Literature Review/Pics/sum_2_Gauss_y_HMC.pdf';
    export_fig(filename,f,'-a1',['-r',num2str(r)]);%, '-grey','-nocrop'
    saveas(f,[filename(1:end-4) '.fig'])
%% Cleaning up working folder
delete anon_model.stan anon_model.hpp
