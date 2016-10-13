function X = state_dynamics( X,theta )
%STATE_DYNAMICS Generating new state particles from the current state
%particles X and the parameter theta.
%
%The specifics of this function depend on the particular dynamical system,
%i.e. this function is user-defined.
%
%This function can be translated to .mex by Matlab Coder. No speed gain
%observed, however. Matlab Coder also does not support bsxfun and .* very
%well. So the code must be explicit in term of dimensionality of variables.

X.state = 0.7*X.state+randn(size(X.state))*sqrt(10^theta);

end

