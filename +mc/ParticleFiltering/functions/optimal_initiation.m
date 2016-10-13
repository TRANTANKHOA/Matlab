function [ X ] = optimal_initiation( N, theta, y, rho )
%OPTIMAL_INITIATION Generate samples from the first posterior distribution
%given the first data point.

Sigma_x = 10^theta/(1-0.7^2);
X = state_initiation(N,theta);
log_prior = @(x) log_normpdf(x,Sigma_x);
log_posterior = @(x) log_normpdf(x,Sigma_x) + log_normpdf(bsxfun(@minus,y,0.5*x),0.1);
X.state = sequential_slice_sampling( X.state, log_prior, log_posterior, rho );

end

