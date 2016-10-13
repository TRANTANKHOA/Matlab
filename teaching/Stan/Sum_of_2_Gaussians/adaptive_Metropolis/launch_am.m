% This script implement the adaptive Metropolis sampler for the model
% 
% $$ \alpha_t \sim N(0,1); t=1,2,...T $$
% 
% $$ y_t \sim N(\alpha_t,\exp(\theta)) $$
% 
%% Initialisation
IN.X = 0;
IN.Sig = 1;            % Initial covariance 
IN.M = M;              % Number of MCMC samples
IN.data = y;           % Necessary data to compute the density function
IN.dsp = 0;            % Display status of MCMC runs, 1 - display | 0 - suppress
IN.cov = 1;            % Adaptive covariance matrix, 1 - adaptive | 0 - fixed
IN.hgrd= 1e2;          % Number of grid point in histograms, must be multiple of 10
IN.filename='am_samples';  % Store MCMC samples in hard drive and save RAM
% target = @(x) logpdf(x,IN.data);
mu_x = 1; sigma_x = 1.5;
target = @(x) log(normpdf(x,mu_x,sigma_x));
%% Adaptive Metropolis
OUT = am(target,IN);
time_am = toc
%% Plotting 
name = '$\upsilon$';
f = figure('name','Adaptive Metropolis');
% trace
subplot(2,2,1)
plot(OUT.TH(:,1))
ylabel(name)
xlabel('Iteration')
axis tight
YL = ylim;
title('Sample trace')
% histograms
subplot(2,2,2)
counts = OUT.hist.p(:,1);
bins = OUT.hist.x(:,1);
barh(bins,counts)
line([0 max(counts)],OUT.X_bar(1)*[1 1],'color','g')
legend('Histogram','Sample mean','Location','southeast')
ylabel(name)
xlabel('Density')
axis tight
set(gca,'Ylim',YL)
title('Sample histogram')
% autocorrelation
subplot(2,2,3)
n_cor = 50;
autocorr(OUT.TH(:,1),n_cor);
hold on
text(mean(xlim),mean(ylim),['$\tau_{\widehat{' name(2:end-1) '}}$ = ' num2str(iact(OUT.TH(:,1)),3)])
title('Sample Correlation')
% zoom in
subplot(2,2,4)
plot(OUT.TH(1:n_cor,1),'Marker','o','LineStyle','-.')
ylabel(name)
xlabel('Iteration')
axis tight
title('Close-up of sample trace')
subplot(2,2,1)
set(gca,'Ylim',YL_theta)
subplot(2,2,2)
set(gca,'Ylim',YL_theta)
drawnow
%% Printing
w = 0.65; h = 0.5; w=round(w*1050); h=round(h*800); r=300; sz = [-1000 1000 w h];
undock(f); set(f,'Units','pixels','Position',sz);
filename = '/Users/Eric/Dropbox/PhD/Reports/Literature Review/Pics/sum_2_Gauss_theta_Metropolis_s.pdf';
export_fig(filename,f,'-a1',['-r',num2str(r)]);%, '-grey','-nocrop'
saveas(f,[filename(1:end-4) '.fig'])
    