function [X0] = state_initiation( N, theta )
%STATE_INITIATION Generate n particles for the state given the parameter
%theta. 
%
%The specifics of this function depend on the particular dynamical system,
%i.e. this function is user-defined.
%
%The unnormalised weights are given on the log-scale as zeros.

%% Generate state particles
sigma_x = 10^theta/(1-0.7^2);
X0.state = sqrt(sigma_x)*randn(N,1);

X0.weight = zeros(N,1);
X0.w = 1/N*ones(N,1);
X0.ESS = 1;
X0.idx = (1:N)';
end

