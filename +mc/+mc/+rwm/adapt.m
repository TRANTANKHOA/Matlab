function target = adapt( target )
%BURNIN Approximate the global covariance matrix and optimal scaling for
%random walk metropolis sampling
%
%   Reference: Andrieu, C., & Thoms, J. (2008). A tutorial on adaptive
%   MCMC. Statistics and Computing.
% _________________________________________________________________________
% INPUT/OUTPUT ARGUMENTS
%
% target        = structure of handles to densities
% target.log    = logarithm of the target density
% target.mu     = mean of the target density
% target.Sigma  = covariance of the target density
% target.chol   = Cholesky factor of Sigma
% target.M      = number of particles required in parallel
% target.Dims   = dimensionality of the target
% target.dirc   = routine for generating auxiliary random direction
% target.map    = routine for mapping to the new proposals
% target.set    = routine for checking set membership of particles
% target.eps    = a scarlar between -0.5 to -1
% target.burnin = initial burnin duration
% target.s      = scaling constant for Metropolis proposal
% target.optacc = theoretically optimal acceptance rate
% target.adp    = adaptation flag
% target.dsp    = display status of burnin runs, 1 - display | 0 - suppress
% target.rep    = desired number of Metropolis sampling repetitions
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

%% Initialise target for RWM sampling
target.adp  = 1;
target      = init_opt(target);
%% Setup burn-in parameters
M           = target.M; target.M = 1;   % Switch to a single chain during burnin
X.pos       = target.mu ;               % Starting location for burnin samples
X.log       = target.log(X.pos);        % Log density of the starting location
X.flags     = 0;                        % Flags for acceptance test
%% Pre-burnin run - generates 'prerun' samples with the initial covariance matrix
for i = 1:target.prerun
    %% Display completion status
    if target.dsp;  target = dsp_status(i,target);  end;
    %% Metropolis sampling
    X       = kernel(target,X);
    %% Update target
    target  = rob_mro(i,target,X);
end
[target.Sigma,target.chol,~] = make_symm(target.Sigma);
%% Burnin run
i=target.prerun+1; target.adp = 0;
while i <= target.rep
    %% Display completion status
    if target.dsp;  target = dsp_status(i,target);  end; 
    %% Metropolis sampling
    X       = kernel(target,X);
    %% Check stopping rule at every chk_stop counts
    if ~mod(i-target.prerun,target.chk_stop)
        if (log(mean(abs(target.delta(i-target.chk_stop+1:i-target.chk_stop/2,1)))) <= log(mean(abs(target.delta(i-target.chk_stop/2:i-1,1))))) &&...
           (log(mean(abs(target.delta(i-target.chk_stop+1:i-target.chk_stop/2,2)))) <= log(mean(abs(target.delta(i-target.chk_stop/2:i-1,2)))))
            target.rep = i;
            target.delta = target.delta(1:target.rep,:);   % Adaptation cost
        end
    end
    %% Update target
    target  = rob_mro(i,target,X);
    i=i+1;
end
target.M = M;
if target.dsp % If requested, give feedback on completion status
    figure;plot(1:target.rep,target.delta(:,2));
    title('Adaptation process of the covariance matrix $$\mathbf{\Sigma}$$');set(gca,'YScale','log');
    xlabel('Number of iterations');ylabel('$$\Delta_m = \Vert\Sigma_m - \Sigma_{m-1}\Vert_1$$');
    axis tight;drawnow
    
    figure;plot(1:target.rep,abs(target.delta(:,1)));
    title('Adaptation process of the scaling factor $$s$$');set(gca,'YScale','log');
    xlabel('Number of iterations');ylabel('$$\Delta_m = \Vert s_m - s_{m-1}\Vert$$');
    axis tight;drawnow
    
    fprintf('Burn-in completed after %d samples\n',target.rep);
    fprintf('Acceptance rate %f %%\n',100*target.acc);
end