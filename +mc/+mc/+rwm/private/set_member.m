function particles = set_member(target,particles)
%SET_MEMBER Checking if the particles are inside the slice or not by
%updating the flags
% _________________________________________________________________________
% INPUT/OUTPUT ARGUMENTS
% 
% target                = logarithm of the target density
% particles.pos         = the positions
% particles.log         = logarithm of the target density at the particles'
%                       positions.
% particles.heit        = the auxiliary random height in slice sampling
% particles.flags       = flags for the vectorisation of Metropolis sampling
% particles.pros        = the position of random proposals
% particles.acc         = acceptance probability of each particle
% particles.log_pros    = log-density of the proposal
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

particles.log_pros  = target(particles.pros);
particles.flags     = particles.log_pros >= particles.heit;
particles.acc       = bsxfun(@min,1,exp(particles.log_pros - particles.log));
particles.log(particles.flags) = particles.log_pros(particles.flags);
