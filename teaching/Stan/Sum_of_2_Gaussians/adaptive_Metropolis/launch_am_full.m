% This script implement the adaptive Metropolis sampler for the model
% 
% $$ \alpha_t \sim N(0,1); t=1,2,...T $$
% 
% $$ y_t \sim N(\alpha_t,\exp(\theta)) $$
% 
%% Initialisation
X = [theta, alpha_T];    % Initial sample
V = eye(T+1);           % Initial covariance
V(1,1) = 1/9;
M = 1e3;                % Number of MCMC samples
OPT.data = y;           % Necessary data to compute the density function
OPT.dsp = 1;            % Display status of MCMC runs, 1 - display | 0 - suppress
OPT.cov = 0;            % Adaptive covariance matrix, 1 - adaptive | 0 - fixed
OPT.hgrd= 1e2;          % Number of grid point in histograms, must be multiple of 10
OPT.filename='am_full_samples';  % Store MCMC samples in hard drive and save RAM
%% Adaptive Metropolis
G = am(X,V,M,OPT);
%% Load back the last batch of MCMC samples to RAM
D = load([G.filename '_' num2str(G.div) '.mat']); % Load TH data
% G.TH = D.TH; 
%% Plotting 
Dims = length(X);
for k=1:2
    if Dims >1
        name = ['$\theta_{' num2str(k) '}$'];
    else
        name = '$\theta$';
    end
    figure
    % trace
    subplot(2,2,1)
    plot(D.TH(:,k))
    line([0 size(D.TH,2)],G.mean(k)*[1 1],'color','g')
    ylabel(name)
    xlabel('Iteration')
    axis tight
    YL = ylim;
    % histograms
    subplot(2,2,2)
    counts = G.hist.p(:,k);
    bins = G.hist.x(:,k);
    barh(bins,counts)
    line([0 max(counts)],G.mean(k)*[1 1],'color','g')
    ylabel(name)
    xlabel('Density')
    axis tight
    set(gca,'Ylim',YL)
    % autocorrelation
    subplot(2,2,3)
    n_cor = 50;
    corr = plot_acf(D.TH(:,k),n_cor);
    hold on
    text(mean(xlim),mean(ylim),['IACT = ' num2str(iact(D.TH(:,k)))])
    % zoom in
    subplot(2,2,4)
    plot(D.TH(k,1:n_cor),'Marker','o','LineStyle','-.')
    ylabel(name)
    xlabel('Iteration')
    axis tight
end
