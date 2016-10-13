function log_p = log_prior_theta(X,theta)
% Calculate the log density of the prior.
SIGMA = 10.^(theta);
log_p = log_normpdf(X(:,:,1),SIGMA/(1-0.7^2),3) + log_normpdf(bsxfun(@minus,X(:,:,2:end),0.7*X(:,:,1:end-1)),SIGMA,3);
end