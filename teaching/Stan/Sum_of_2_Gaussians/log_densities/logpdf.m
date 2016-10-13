function log_dens = logpdf(X,Y)
%% posterior of theta
sigma_2     = 1 + exp(X);
T           = length(Y);
log_dens    = log_prior(X) - 0.5.*( sum(Y.^2./sigma_2) + T.*log(sigma_2) );

%% posterior of (theta,alpha)
% tau_2 = exp(X(1));
% alpha_T = X(2:end);
% T           = length(Y);
% log_dens    = log_prior(X(1)) - 0.5.*( sum(alpha_T.^2 + (Y-alpha_T).^2./tau_2) + T.*log(tau_2) );

end