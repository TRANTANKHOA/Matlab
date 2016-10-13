function [X] = state_propagation( X, y, rho )
%STATE_PROPAGATION Updating the weights by the likelihood function
%
%The specifics of this function depend on the particular dynamical system,
%i.e. this function is user-defined.

%% Calculate the weights by likelihood function
X.weight = log_normpdf(y-0.5*X.state,0.1) + X.weight;

%% Resampling if necessary
[X.w,X.ESS] = weight_normalisation( X.weight );
N=length(X.weight);
if X.ESS<=rho
    [X.state, X.idx] = important_resampling( X.state,X.w,N );
    X.weight(:) = 0;
    X.w(:) = 1/N;
%     X.ESS = 1;
else
    X.idx = (1:N)';
end
end

