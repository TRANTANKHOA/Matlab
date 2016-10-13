function [w, ESS] = weight_normalisation( log_omega )
%WEIGHT_NORMALISATION Normalising the weights and converting them from the
%log-scale to linear-scale. 
%
% This function gives an error when the effective sample size is zero.

% Calculating the effective sample size
sum_omega = lse(log_omega);
if sum_omega ~= -Inf   % Checking of the sum of weights is not 0.
    log_omega = log_omega-sum_omega;   % Normalising the weights
    w = exp(log_omega); w = w/sum(w);  
    ESS = exp(-lse(log_omega*2))/length(w);            % ESS = 1/sum(w.^2)
else
    error('Zero effective sample size');         % Zero effective sample size
end

end

