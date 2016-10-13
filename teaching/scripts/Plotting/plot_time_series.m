function h = plot_time_series(theta,name,n_cor,window)
%% Plotting theta
Dims = size(theta,2);
for k=1:Dims
    h(k) = figure('name',window);
    % trace
    subplot(2,2,1)
    plot(theta(:,k))
    ylabel(name)
    xlabel('Iteration')
    axis tight
    YL = ylim;
    title('Sample trace')
    % histograms
    subplot(2,2,2)
    [counts,bins] = hist(theta,50); 
    counts= counts/(sum(counts)*(bins(2)-bins(1)));
    barh(bins,counts)
    line([0 max(counts)],mean(theta(:,k))*[1 1],'color','g')
    legend('Histogram','Sample mean','Location','northeast')
    ylabel(name)
    xlabel('Density')
    axis tight
    set(gca,'Ylim',YL)
    title('Sample histogram')
    % autocorrelation
    subplot(2,2,3)
    autocorr(theta(:,k),n_cor);
    hold on
    text(mean(xlim),mean(ylim),['$\tau_{\widehat{' name(2:end-1) '}}$ = ' num2str(iact(theta(:,k)),3)])
    title('Sample Correlation')
    % zoom in
    subplot(2,2,4)
    plot(theta(1:n_cor,k),'Marker','o','LineStyle','-.')
    ylabel(name)
    xlabel('Iteration')
    axis tight
    title('Close-up of sample trace')
end