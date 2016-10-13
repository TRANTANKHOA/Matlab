% This script generates the posterior samples for the model
% 
% $$ \alpha_t \sim N(\alpha_{t-1}\exp(-\phi\alpha_{t-1}^2),\exp(\theta)); t=1,2,...T $$
% 
% $$ y_t \sim N(\alpha_t,1) $$
% 
%% Defining true parameters
clc; close all
theta_0 = 0; phi_0 = 2; tau_0 = sqrt(exp(theta_0));
%% Creating simulated data
T       = 200;                  % data points
alpha   = randn(1,T)*tau_0;
for t=2:T
    alpha(t) = alpha(t) + exp(-phi_0*alpha(t-1)^2)*alpha(t-1);
end
y       = phi_0*exp(alpha) + randn(1,T);
plot([y' alpha'])
%% Sampling $\theta$ using Hamiltonian Monte Carlo
tic
launch_stan
time_stan = toc
%% Cleaning up working folder
delete *.stan *.hpp *.mat *.csv *.R


