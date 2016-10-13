function [ Path ] = auxiliary_particle_filtering(  N, theta, y, rho )
%AUXILIARY_PARTICLE_FILTERING Perform auxiliary particle filtering based on 
%sequential slice sampling.
%
%INPUTS:
% N     - the number of desired particles
% theta - model parameter
% y     - record of output data
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
%% Sample from the first posterior
Path = optimal_initiation( N, theta, y(1), rho );
Path = repmat(Path,T,1);
for i = 2:T
    Path(i) = optimal_proposal( Path(i-1), N, theta, y(i), rho );
end

