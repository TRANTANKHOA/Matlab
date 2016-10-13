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
    window = 'Hamiltonian Monte Carlo';
    f = plot_time_series(fit.extract('permuted',true).theta,'$\upsilon$',10,window);
    subplot(2,2,1)
    YL_theta = ylim;
    subplot(2,2,2)
    legend('off')

    %% Printing
    w = 0.5; h = 0.35; w=round(w*1050); h=round(h*800); r=300; sz = [-1000 1000 w h];
    drawnow; undock(f); set(f,'Units','pixels','Position',sz);
    filename = '/Users/Eric/Dropbox/PhD/Reports/Literature Review/Pics/sum_2_Gauss_theta_HMC.pdf';
    export_fig(filename,f,'-a1',['-r',num2str(r)]);%, '-grey','-nocrop'
    saveas(f,[filename(1:end-4) '.fig'])

    %% Plotting generated data y
    Y = fit.extract('permuted',true).y_tilde;
    plot_y_T
    YL_y = ylim;
    %% Printing
    w = 0.35; h = 0.5; w=round(w*1050); h=round(h*800); r=300; sz = [-1000 1000 w h];
    f = gcf; undock(f); set(f,'Units','pixels','Position',sz);drawnow; 
    filename = '/Users/Eric/Dropbox/PhD/Reports/Literature Review/Pics/sum_2_Gauss_y_flat.pdf';
    export_fig(filename,f,'-a1',['-r',num2str(r)]);%, '-grey','-nocrop'
    saveas(f,[filename(1:end-4) '.fig'])
    
    %% Cleaning up working folder
    delete anon_model.stan anon_model.hpp
