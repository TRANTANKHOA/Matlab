function particles = bisection(particles)
%BISECTION Tunning the fictional tempering parameter lambda using bisection 
% _________________________________________________________________________
% INPUT/OUTPUT ARGUMENTS
%
% particles             = the particles
% particles.pos         = the positions
% particles.w           = the weights
% particles.particles.log_w       = logarithm of the weights
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

new_lbda = 1;   % Start with the most ambitious increment for lambda
%% Bisection
while(true)  
    %% Calculating the next density values
    new_trst = lwse_row([particles.log particles.inst],[new_lbda (1-new_lbda)]);
    %% Calculating the weights on the log scale
    particles.log_w = new_trst-particles.trst;
    %% Calculating the log of effective sample size
    log_sum = lse(particles.log_w);
    if log_sum ~= -Inf   % Checking of the sum of weights is not 0.
        particles.log_w = particles.log_w-log_sum;  % Normalising the weights
        log_ESS = -lse(particles.log_w*2);% ESS = 1/sum(w.^2)
    else
        log_ESS = -Inf;         % Zero effective sample size
    end
    %% Evaluate the condition for stoping Bisection
    if log_ESS >= particles.threshold
        break;  % Escaping the loop with a new value for lambda
    else
        if (new_lbda-particles.lambda) < 1e-6
%             keyboard
            error('Bisection failed');      % The problem is too hard for automatic Monte Carlo sampling
        else
            new_lbda = (new_lbda+particles.lambda)/2;   % Repeating bisection with new lambda
        end
    end
end
%% Normalising the importance weights
[particles.w, particles.ESS] = mc.is.weight_normalisation( particles.log_w );
%% Updating output
particles.trst    = new_trst;
particles.lambda  = new_lbda;