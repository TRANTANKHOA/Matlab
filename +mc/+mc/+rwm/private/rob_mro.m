function target = rob_mro(idx,target,particles)
%rob_mro Tuning the scale and covariance of random walk proposals using
%Robbins-Monro recursion
%
% _________________________________________________________________________
% INPUT/OUTPUT ARGUMENTS
%
% idx           = current index
% target.mu     = mean of the target density
% target.Sigma  = covariance of the target density
% target.chol   = Cholesky factor of Sigma
% target.Dims   = dimensionality of the target
% target.s      = scaling constant for Metropolis proposal
% target.optacc = theoretically optimal acceptance rate
% target.eps    = a scarlar between -0.5 to -1
% target.dirc   = routine for generating auxiliary random direction
% target.map    = routine for mapping to the new proposals
% target.set    = routine for checking set membership of particles
% target.delta  = step change in adaptation
% target.acc    = average acceptance rate
% target.pct    = percentage complete count initialised to zero;
% target.mark   = used to keep track of elapsed time
% target.adp    = adaptation flag
% target.chk_stop   = periodically check to stop burn-in
% target.prerun     = number of samples before burnin
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
%   Created:        07-Aug-2016
%   Last edited:    07-Aug-2016
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

%% Update scale
target.mu   = target.mu*(idx-1)/idx + particles.pos/idx;
target.acc  = target.acc*(idx-1)/idx + particles.acc/idx;
gamma       = idx^target.eps;
target.delta(idx,1)  = gamma*(particles.acc-target.optacc);
target.s = sqrt(exp(2*log(target.s)+target.delta(idx,1)));
%% Update Sigma
delta_X = particles.pos-target.mu;
delta_S = (delta_X'*delta_X - target.Sigma)./idx;
target.delta(idx,2) = norm(delta_S,1);
target.Sigma = target.Sigma+delta_S;
if target.adp == 0
    [target.Sigma,target.chol,~] = make_symm(target.Sigma);
end