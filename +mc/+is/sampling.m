function particles = sampling(target)
%IMPORTANCE_SAMPLING generate M weighted samples to approximate the target
% density by sampling from the instrm density function, where as:
% *target* is the function handle for the log density function
% *target.inst* is the handle to the importance density
% _________________________________________________________________________
% INPUT ARGUMENTS
%
% target        = structure of handles to densities
% target.log    = logarithm of the target density
% target.inst   = logarithm of the instrumental density
% target.rand   = random number generator for the instrumental density
% target.M      = number of particles required
% _________________________________________________________________________
% OUTPUT ARGUMENTS
% 
% particles         = the particles
% particles.pos     = the positions
% particles.w       = the weights
% particles.log_w   = logarithm of the weights
% particles.ESS     = effective sample sizes
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

%% Generate random vectors from the instrumental density
particles.pos = target.rand(target);

%% Calculating the weights on the log scale
particles.log_w = target.log(particles.pos) - target.inst(particles.pos);

%% Normalising the weights into the linear scale and compute the effective sample size.
[particles.w,particles.ESS] = mc.is.weight_normalisation( particles.log_w );
