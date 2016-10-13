%% This script generates the posterior samples for the model
% 
% $$ x_t \sim N(\phi x_{t-1},\tau^2); t=1,2,...T $$
% 
% $$ y_t \sim N(x_t,\sigma^2) $$
% 
%% Defining true parameters
clc; close all;clear
sys.sigma   = .1;
sys.tau     = 1; 
sys.phi     = 0.9; 
sys.k       = sys.tau/sqrt(1-sys.phi^2);
%% Creating simulated data
sys.T           = 200;  % data points
sys.x(sys.T,1)  = 0;    % allocate memory
sys.x(1)        = randn*sys.k;
for t=2:sys.T
    sys.x(t) = sys.phi*sys.x(t-1) + randn*sys.tau;
end
sys.y           = sys.x + randn(sys.T,1)*sys.sigma;  % Y = X + e
%% Plot simulated data
figure(1);
[hAx,hLine1,hLine2] = plotyy(1:sys.T,sys.x,1:sys.T,sys.y,'plot','stem');
title('Simulated Output and True State Trajectories')
xlabel('$t$')
ylabel(hAx(1),'$x_t$') % left y-axis
ylabel(hAx(2),'$y_t$') % right y-axis
hLine1.LineStyle    = '-.';
hLine1.Marker       = '*';
set(hAx(1),'Ylim',get(hAx(2),'Ylim'))
%% Matrix representation of linear state space model
sys.A       = eye(sys.T); sys.A(2:sys.T+1:sys.T*(sys.T-1)) = -sys.phi;
sys.g       = zeros(sys.T,1); sys.g(1) = 1;
sys.Va      = sys.tau^2*eye(sys.T); sys.Va(1,1) = 0;
sys.Ainv_g  = sys.A\sys.g;
sys.Vx      = sys.Ainv_g*sys.Ainv_g'*sys.k^2 + sys.Va;
sys.Vy      = sys.Vx + sys.sigma^2*eye(sys.T);
Exy     = sys.Vx*((sys.Vy)\sys.y);
Vxy     = sys.Vx - sys.Vx*(sys.Vy\sys.Vx);
xtt     = Exy;
stt     = diag(Vxy);
%%
figure(2);
[hAx,hLine1,hLine2] = plotyy(1:sys.T,sys.x,1:sys.T,xtt);
title('Estimated and True State Trajectories')
xlabel('$t$')
ylabel(hAx(1),'$x_t$') % left y-axis
ylabel(hAx(2),'$E(x_t \mid y_{1:t})$') % right y-axis
hLine1.LineStyle    = '-.';
hLine1.Marker       = '*';
hLine2.LineStyle    = '-';
hLine2.Marker       = 'o';
set(hAx(2),'Ylim',get(hAx(1),'Ylim'))
axis(hAx(2));
hold on
plot(1:sys.T,xtt+3*sqrt(stt),'-y')
plot(1:sys.T,xtt-3*sqrt(stt),'-y')
%% Maximum likelihood estimation for sigma
n       = 200; 
theta   = linspace(-6,6,n)'; % \theta = \log(\sigma^2)
l(n,1)  = 0; % memory allocation
func    = @(theta) likelihood(sys,[theta log(sys.tau^2)]);
for i = 1:n
    l(i) = func(theta(i));
end
est_theta = fminunc(func,0);
est_sigma = sqrt(exp(est_theta));
figure(3)
plot(theta,l);xlabel('Parameter - $\theta = \log(\sigma^2)$');ylabel('Negative of Log-likelihood');
line(est_theta*[1 1],ylim,'LineStyle','-.')
text(mean(xlim),mean(ylim),['$\hat{\sigma} = ' num2str(est_sigma,4) '$'])
%% Maximum likelihood estimation for tau
n       = 200; 
theta   = linspace(-6,6,n)'; % \theta = \log(\tau^2)
l(n,1)  = 0;
func    = @(theta) likelihood(sys,[log(sys.sigma^2) theta]);
for i = 1:n
    l(i) = func(theta(i));
end
est_theta = fminunc(func,0);
est_tau   = sqrt(exp(est_theta));
figure(4)
plot(theta,l);xlabel('Parameter - $\theta = \log(\tau^2)$');ylabel('Negative of Log-likelihood');
line(est_theta*[1 1],ylim,'LineStyle','-.')
text(mean(xlim),mean(ylim),['$\hat{\tau} = ' num2str(est_tau,4) '$'])

%% Maximum likelihood estimation for (sigma,tau)
func    = @(theta) likelihood(sys,theta);
est_theta   = fminunc(func,2*[sys.sigma sys.tau]);
n       = 30; w = 2;
x       = linspace(est_theta(1)-w,est_theta(1)+w,n)'; % x = \log(\sigma^2)
y       = linspace(est_theta(2)-w,est_theta(2)+w,n)'; % y = \log(\tau^2)
figure(5);plot3D(x,y,func);
est_sigma   = sqrt(exp(est_theta(1)));
est_tau     = sqrt(exp(est_theta(2)));
line(est_theta(1)*[1 1],est_theta(2)*[1 1],zlim,'LineStyle','-.')
xlabel('Parameter - $\theta_1 = \log(\sigma^2)$');
ylabel('Parameter - $\theta_2 = \log(\tau^2)$');
zlabel('Negative of Log-likelihood');
%% Sampling using Hamiltonian Monte Carlo
tic
launch_stan
time_stan = toc
%% Cleaning up working folder
delete *.stan *.hpp *.mat *.csv *.R anon_model 


