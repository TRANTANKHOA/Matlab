function log_dens = log_prior(theta) % log of prior density of the model
sig_tht     = 1;    
mu_tht      = 0;
log_dens    = log_normpdf(theta-mu_tht,sig_tht); 
end
