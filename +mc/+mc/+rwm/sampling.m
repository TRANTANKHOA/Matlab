function particles = sampling( target,particles )
%RWM_SAMPLING Template for Gaussian random walk sampling algorithm. 
%
% _________________________________________________________________________
% INPUT ARGUMENTS
%
% target        = structure of handles to densities
% target.log    = logarithm of the target density
% target.mu     = mean of the target density
% target.Sigma  = covariance of the target density
% target.chol   = Cholesky factor of Sigma
% target.inst   = logarithm of the instrumental density
% target.rand   = random number generator for the instrumental density
% target.M      = number of particles required in parallel
% target.Dims   = dimensionality of the target
% target.name   = name for saving file to hard disk
% target.dirc   = routine for generating auxiliary random direction
% target.map    = routine for mapping to the new proposals
% target.set    = routine for checking set membership of particles
% target.m      = number of repetitions per division
% target.div    = number of division
% target.rep    = number of repetitions
% target.eps    = a scarlar between -0.5 to -1
% target.burnin = initial burnin duration
% target.s      = scaling constant for Metropolis proposal
% target.optacc = theoretically optimal acceptance rate
% target.mod    = mode for generating auxiliary innovation
% _________________________________________________________________________
% OUTPUT ARGUMENTS
%
% particles         = the particles
% particles.pos     = the positions
% particles.log     = logarithm of the target density at the particles'
%                     positions.
% particles.heit    = the auxiliary random height 
% particles.dirc    = the auxiliary random innovation 
% particles.flags   = flags for the accept/reject status
% particles.pros    = the position of random proposals
%
% -------------------- Copyright (C) 2016 Khoa T. Tran --------------------
%   Created:        04-Aug-2016
%   Last edited:    04-Aug-2016
%   Matlab version: 8.6.0.267246 (R2015b)
%	Email:          khoa.tran@unsw.edu.au
%
%   This program is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation, either version 3 of the License, or (at
%   your option) any later version.
%       This program is distributed in the hope that it will be useful, but
%   WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
%   General Public License for more details.
%       A copy of the GNU General Public License can be found at:
%   http://www.gnu.org/licenses/.
% _________________________________________________________________________

%% Initialise target for RWM sampling
target.adp  = 0;
target      = init_opt(target);
%% Compute the Cholesky factor if none is found
if ~isfield(target, 'chol') || isempty(target.chol)
    [target.Sigma,target.chol,~] = make_symm(target.Sigma);
end
%% Initialise particles for RWM sampling
particles.pos = target.rand(target);
particles.log = target.log(particles.pos);
if target.mod 
    target.pros     = @(v)log_mvnpdf(v,zeros(1,target.Dims),target.s^2*target.Sigma);   
    particles.dirc  = zeros(target.M,target.Dims);
    particles.auxlog= target.pros(particles.dirc);
    target.kernel   = @s_kernel;
else
    particles       = kernel(target,particles);
    target.kernel   = @kernel;
end
%% Sampling
for count = 1:target.div
    for k = 1:target.m
        particles(1+k) = target.kernel(target,particles(k)); 
    end
    %% Save SMC samples to hard disk for retrospection
    particles(1) = [];
    save(['random_walk_' target.name '_' num2str(count)],'particles'); 
    particles(1) = particles(end);
end