function particles = sampling( target )
%SMC Sequential Monte Carlo algorithm for static target density. 
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
% 
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


%% Generate N samples from the instrumental density
transition          = target;    
particles.pos       = target.rand(target);        
particles.lambda    = 0; 
particles.ESS       = 0;
particles.threshold = log(0.5*target.M);            % Setting the threshold for bisection algorithm
count               = 1;
while particles.lambda<1
    %% Bisection
    particles.inst  = target.inst(particles.pos);   % Storing the values of instrumental density 
    particles.log   = target.log(particles.pos);    % Storing the values of target density
    particles.trst  = lwse_row([particles.log particles.inst],[particles.lambda (1-particles.lambda)]); % Current density
    particles       = bisection(particles);
    %% Messaging   
    clc
    formatSpec = 'Performing Sequential Slice Sampling.\nEffective sample size is %2.1f %% samples.\n';
    fprintf(formatSpec,particles.ESS/target.M*100)
    formatSpec = 'Lambda = %1.6f.\n';
    fprintf(formatSpec,particles.lambda)
%     keyboard
    %% Importance Resampling
    idx = mc.is.resampling( particles.w, target.M );
    % Retabulate the values of particles and its density
    particles.pos   = particles.pos(idx,:);
    particles.log   = particles.trst(idx,:);
    %% Executing an MCMC Kernel, e.g. slice sampler 
    % Define the current target 
    transition.log  = @(x) lwse_row([target.log(x) target.inst(x)],[particles.lambda (1-particles.lambda)]);
    transition.dirc = @(particles) mc.mc.slice.dirc.rnd_dir(particles.pos);
    transition.map  = @(theta,v,r) theta+bsxfun(@times,v,r); % x = theta+v.r
    transition.set  = @(particles) mc.mc.slice.set_member(transition.log,particles);
    % Recursive sampling
    particles = mc.mc.slice.recrsv(transition, particles);
    %% Save SMC samples to hard disk for retrospection
    save(['smc_' target.name '_' num2str(count)],'particles','-v7.3'); count = count + 1;
end

end

