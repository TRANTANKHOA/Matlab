function target = init_opt(target)
%init_opt Setting up initial parameters for adaptive Metropolis according
% to Roberts, G. O., & Rosenthal, J. S. (2001). Optimal Scaling for Various
% Metropolis-Hastings Algorithms. Statistical Science, 16(4), 351?367.
% Retrieved from http://projecteuclid.org/euclid.ss/1015346320
% _________________________________________________________________________
% INPUT/OUTPUT ARGUMENTS
%
% target.Dims   = dimensionality of the target
% target.s      = scaling constant for Metropolis proposal
% target.optacc = theoretically optimal acceptance rate
% target.eps    = a scarlar between -0.5 to -1
% target.mod    = mode for generating auxiliary innovation
% target.dirc   = routine for generating auxiliary random direction
% target.map    = routine for mapping to the new proposals
% target.set    = routine for checking set membership of particles
% target.delta  = step change in adaptation
% target.acc    = average acceptance rate
% target.pct    = percentage complete count initialised to zero;
% target.mark   = used to keep track of elapsed time
% target.chk_stop = periodically check to stop burn-in
% target.prerun   = number of samples before burnin
% -------------------- Copyright (C) 2016 Khoa T. Tran --------------------
%   Created:        05-Aug-2016
%   Last edited:    05-Aug-2016
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

if target.adp 
    switch target.Dims              % Define optimal acceptance rate
        case 1
            target.optacc=0.441;
        case 2
            target.optacc=0.352;
        case 3
            target.optacc=0.316;
        case 4
            target.optacc=0.285;
        case 5
            target.optacc=0.275;
        case 6
            target.optacc=0.273;
        case 7
            target.optacc=0.270;
        case 8
            target.optacc=0.268;
        case 9
            target.optacc=0.267;
        case 10
            target.optacc=0.266;
        case 11
            target.optacc=0.265;
        case 12
            target.optacc=0.264;
        otherwise
            target.optacc=0.234;
    end
    target.s    = 2.38/sqrt(target.Dims);   % Set initial scaling factor
    target.eps  = -0.8;                     % Could be -0.5 to -1
    target.delta= zeros(target.rep,2);      % Step change in adaptation
    target.acc  = 0;                        % Average acceptance rate
    target.pct  = 0;                        % Percentage complete count initialised to zero;
    target.mark = tic;              % Used to keep track of elapsed time
    target.chk_stop = 2000;         % Periodically check to stop burn-in
    target.prerun   = ceil(target.rep/1e2); % Some initial samples are required before the empirical covariance stabilised
end
%% Define the current target 
% Subroutine for getting auxiliary innovation vectors
target.dirc = @(target) target.s .* randn(target.M,target.Dims) * target.chol;
% Subroutine for mapping the new proposals
target.map  = @(particles) particles.pos+particles.dirc; % x = theta+v
% Subroutine for checking set membership
target.set  = @(particles) set_member(target.log,particles);
