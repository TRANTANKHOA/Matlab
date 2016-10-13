function [ Path ] = particle_filtering(  N, theta, y, rho )
%PARTICLE_FILTERING Perform particle filtering based on sequential
%important sampling and resampling.
%
%INPUTS:
% N     - the number of desired particles
% y     - record of output data
% theta - model parameter
% rho   - resampling threshold
%
%OUTPUTS:
% Path        - Array of filtering distribution stored in a structure with
% Path.state  - State particles from the sequence of filtering distributions
% Path.weight - Particle weights on the log-scale
%
%This function can be translated to .mex by Matlab Coder. No speed gain
%observed, however. Matlab Coder also does not support bsxfun and .* very
%well. So the code must be explicit in term of dimensionality of variables.

T = length(y);      

Path = state_propagation(state_initiation(N,theta), y(1), rho);
Path = repmat(Path,T,1);

for i = 2:T
    Path(i) = state_propagation(state_dynamics(Path(i-1),theta), y(i), rho);
end

