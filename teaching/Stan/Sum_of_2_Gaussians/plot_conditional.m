% This script describe the behaviour of the conditional density of theta
% given alpha 
%% Setup range of plot
clearvars; close all;clc
theta = (-2:0.01:0)';
tau_2 = exp(theta);
%% Generate experiment data
T = 100;
alpha = rand(1,T);
y = alpha + rand(1,T);
%% Plotting the conditional density of theta
figure
subplot(1,2,1)
hold on
xlabel('$\theta$')
ylabel('$p(\theta | \alpha, Y_T)$')
subplot(1,2,2)
hold on
xlabel('$\tau^2$')
ylabel('$p(\tau^2 | \alpha, Y_T)$')

for k=1:1000
alpha = rand(1,T); % randomise alpha
% y = alpha + rand(1,T);
log_density = -0.5.*sum( bsxfun(@plus,bsxfun(@rdivide,(y-alpha).^2,tau_2), alpha.^2),2 )...
- log(tau_2.^(T/2)) + log_prior(theta) ;
subplot(1,2,1)
plot(theta,(exp(log_density)))
subplot(1,2,2)
% p(tau) = p(theta)d.theta/d.tau
plot(tau_2,exp(log_density)./tau_2)
end