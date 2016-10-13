function log_l = log_likelihood(X,y)
% Calculate the log density of the prior.

log_l = log_normpdf(bsxfun(@minus,y,0.5*X(:,2:end)), 0.1 );

end