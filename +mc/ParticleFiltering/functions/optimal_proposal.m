function [ X1 ] = optimal_proposal( X0, N, theta, y, rho )
%OPTIMAL_PROPOSAL Generate samples from the optimal proposal given the
%previous state estimation
Xm = mean(X0.state,1);

Sigma_x = 10^theta;
X1.state = 0.7*Xm + randn(N,1)*sqrt(Sigma_x);
log_prior = @(x) log_normpdf(x-0.7*Xm,Sigma_x);
log_posterior = @(x) log_normpdf(x-0.7*Xm,Sigma_x) + log_normpdf(bsxfun(@minus,y,0.5*x),0.1);
X1.state = sequential_slice_sampling( X1.state, log_prior, log_posterior, rho );
X1.weight = X0.weight + log_normpdf(X1.state-0.7*X0.state,Sigma_x) - log_prior(X1.state);

%% Resampling if necessary
[X1.w,X1.ESS] = weight_normalisation( X1.weight );
if X1.ESS<=rho
    [X1.state, X1.idx] = important_resampling( X1.state,X1.w,N );
    X1.weight(:) = 0;
    X1.w(:) = 1/N;
%     X.ESS = 1
else
    X1.idx = (1:N)';
end
end