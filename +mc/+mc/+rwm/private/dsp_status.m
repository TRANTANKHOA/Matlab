function target = dsp_status(idx,target)
%dsp_status Display completion status for MCMC run
%
% _________________________________________________________________________
% INPUT/OUTPUT ARGUMENTS
% idx         = current index
% target.rep  = number of iterations
% target.pct  = percentage of completion
% target.mark = time target.marker
%
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
if mod(idx,floor(target.rep/20))==0
    target.pct=target.pct+5;
    elaps = toc(target.mark);
    fprintf('Complete = %d%% burn-in, Time since last update = %f s\n',target.pct,elaps); 
    if idx>2
    remaining = (20 - idx/floor(target.rep/20))*(elaps);
    hrs = floor(remaining/3600); remaining = rem(remaining,3600);
    mins = floor(remaining/60); 
    secs = floor(rem(remaining,60));
    fprintf('Predicted burn-in completion in %d hrs:%d mins:%d secs\n',hrs,mins,secs)
    fprintf('Acceptance rate %f %%\n',100*target.acc);
    end;
    target.mark=tic;
end;