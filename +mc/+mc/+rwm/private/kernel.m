function particles = kernel(target,particles)
%METROPOLIS_KERNEL One step Metropolis sampling. 
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
% particles.acc     = acceptance probability of each particle
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

%% Generate random innovation vectors
particles.dirc  = target.dirc(target); 
%% Draw a random height on the log scale for each sample.
particles.heit  = particles.log-exprnd(1,target.M,1);  
%% Mapping out the proposals
particles.pros  = target.map(particles);
%% Accept - reject test
particles       = target.set(particles);    
%% Update samples
particles.pos(particles.flags,:) = particles.pros(particles.flags,:);
