function log_p = log_prior(X)
% Calculate the log density of the prior.

log_p = log_normpdf(X(:,1)) + log_normpdf(bsxfun(@minus,X(:,3:end),0.7*X(:,2:end-1)), 10.^X(:,1) );

end