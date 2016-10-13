% This script generates the posterior samples for the model
% 
% $$ \alpha_t \sim N(0,1); t=1,2,...T $$
% 
% $$ y_t \sim N(\alpha_t,\exp(\theta)) $$
% 
%% Creating simulated data
clc; close all; clearvars 
rng default
T       = 100;                  % data points
alpha   = randn(1,T);
y       = alpha + randn(1,T);
%% Setting up initial samples
theta   = 0;    % first parameter
alpha_T = y./2; % mean value of the Kalman filter output
%% Sampling $\theta$ using Hamiltonian Monte Carlo
tic
launch_stan
time_stan = toc
stan_exit;
drawnow
%% Sampling $\theta$ using adaptive Metropolis sampler
tic
launch_am
%% Sampling $(\theta,\alpha_t)$ using T-factors Auxiliary Gibbs sampler
tic
Ideal_Gibbs
%% Sampling $(\theta,\alpha_t)$ using Hamiltonian Monte Carlo
tic
launch_stan_full
time_stan_full = toc
stan_exit_full;
drawnow
%% Sampling $(\theta,\alpha_t)$ using T-factors Auxiliary Gibbs sampler
tic
Auxiliary_Gibbs
%% Sampling $(\theta,\alpha_t)$ using slice sampler
tic
my_slice
%% Plotting the cpu time
time = [time_stan_full time_I_Gibbs time_slice time_Gibbs];
names = {'HMC-C++'; 'Gibbs-Matlab'; 'Slice-Matlab'; 'Aux. Gibbs-Matlab'};
f = figure('name','Computing Time');
bar(time)
% set(gca, 'XTickLabel',names, 'XTick',1:numel(names))
set(gca,'xticklabel',names)
xlabel('Algorithms');ylabel('Computing Time')
    %% Printing
%     set(gca,'Ylim',[0,150])
    w = 0.5; h = 0.25; w=round(w*1050); h=round(h*800); r=300; sz = [-1000 1000 w h];
    drawnow; undock(f); set(f,'Units','pixels','Position',sz);
    filename = '/Users/Eric/Dropbox/PhD/Reports/Literature Review/Pics/sum_2_Gauss_CPU.pdf';
    export_fig(filename,f,'-a1',['-r',num2str(r)],'-nocrop');%, '-grey'
    saveas(f,[filename(1:end-4) '.fig'])
%% Cleaning up working folder
delete *.stan *.hpp *.mat *.csv *.R


