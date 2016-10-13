function log_l = log_likelihood_theta(X,y)
% Calculate the log density of the prior.

log_l = log_normpdf(bsxfun(@minus,y,0.5*X(:,:,1:end)), 0.1 ,3);

end