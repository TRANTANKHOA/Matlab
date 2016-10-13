function particles = recrsv(target,particles) 
%RECRSV Template for recursive sampling (i.e. slice sampling kernel). 
% _________________________________________________________________________
% INPUT ARGUMENTS
%
% target        = structure of handles to densities
% target.log    = logarithm of the target density
% target.mu     = mean of the target density
% target.Sigma  = covariance of the target density
% target.inst   = logarithm of the instrumental density
% target.rand   = random number generator for the instrumental density
% target.M      = number of particles required
% target.Dims   = dimensionality of the target
% target.name   = name for saving file to hard disk
% target.dirc   = routine for generating auxiliary random direction
% target.map    = routine for mapping to the new proposals
% target.set    = routine for checking set membership of particles
% _________________________________________________________________________
% OUTPUT ARGUMENTS
% 
% particles             = the particles
% particles.pos         = the positions
% particles.w           = the weights
% particles.log_w       = logarithm of the weights
% particles.ESS         = effective sample sizes
% particles.log         = logarithm of the target density at the particles'
%                       positions.
% particles.inst        = logarithm of the instrumental density at the
%                       particles' positions.
% particles.trst        = logarithm of the transitional density at the
%                       particles' positions.
% particles.lambda      = fictional tempering parameter
% particles.threshold   = lower bound of the effective sample size
% particles.heit        = the auxiliary random height in slice sampling
% particles.dirc        = the auxiliary random direction in slice sampling
% particles.b/a         = the boundaries of the unidimensional slice
% particles.count       = counting the number of density evaluation per
%                       slice
% particles.flags       = flags for the vectorisation of parallel slice
%                       sampling
% particles.pros        = the position of random proposals
%
% -------------------- Copyright (C) 2016 Khoa T. Tran --------------------
%   Created:        31-Jul-2016
%   Last edited:    31-Jul-2016
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


% Generate random vector v for each sample
particles.dirc  = target.dirc(particles); 
% Draw a random height on the log scale for each sample.
particles.heit  = particles.log-exprnd(1,target.M,1);   
% Creating first interval
expand;
% Proposing and shrinking bounds
shrink;
end
