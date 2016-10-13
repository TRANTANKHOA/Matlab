function idx = resampling( w,N )
%IMPORTANCE_RESAMPLING Perform a systematic resampling and return the index
%
% _________________________________________________________________________
% INPUT ARGUMENTS
%
% w     = the weights
% N     = number of particles required
% _________________________________________________________________________
% OUTPUT ARGUMENTS
%
% idx     = the indexes
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

u = rand/N;           % Draw u uniformly from [0,1/ESS]
S = w(1);             % Set S = omega_1              
idx = zeros(N,1);     % Create storage for the indexes

% Sequentially filling the value of the idx vector
k = 1;                  
for i=1:N
    while S <= u
        k = k+1;
        S = S+w(k);   % Increase S
    end
    idx(i) = k;       % Store the new index 
    u = u+1/N;        % Increase u
end
