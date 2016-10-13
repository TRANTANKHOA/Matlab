function [target,particles] = adaptive_sampling( target )
%ADAPTIVE_IMPORTANCE_SAMPLING Tuning the instrumental density and perform
%importance sampling
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
%	  Email:          khoa.tran@unsw.edu.au
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

%% Calling importance sampling
particles   = mc.is.sampling(target);
%% Initialise a Do-While loop
s = 0;  new_target = target;
while(true)                
    %% Updating the mean value
    new_target.mu = sum(bsxfun(@times,particles.pos,particles.w));
    %% Updating the covariance only when ESS > 3
    if particles.ESS > 3*target.Dims
        new_target.Sigma        = weightedcov(particles.pos,particles.w);
        new_target.M            = target.M+round(particles.ESS);
        %% Calculating the Cholesky factor of Sigma to check positive definiteness
        [new_target.Sigma,~,~]  = make_symm(new_target.Sigma);
    else
        new_target.Sigma        = 10*eye(target.Dims)*0.9^s;
    end
    %% Call importance sampling to generate new weighted samples
    new_target.inst = @(x) log_mvnpdf(x,new_target.mu,new_target.Sigma);
    new_particles   = mc.is.sampling( new_target );
    %% Evaluate the termination of the loop
    if (new_particles.ESS/new_target.M < particles.ESS/target.M)  % new result is worse
        if (particles.ESS > 3*target.Dims)
            break;  % end the tuning loop
        else
            s = s+1;
            if s > 1e2
                error('Adaptive Importance Sampling failed')
            end
        end
    else            % Storing the last update
        particles   = new_particles;
        target      = new_target;
    end  
    
    if particles.ESS/target.M > 0.5
        break
    end   
end