function v = rnd_dir(particles)
%RND_DIR Get a number of particles and return the difference between two
%randomly selected particles.
% _________________________________________________________________________
% INPUT ARGUMENTS
% 
% particles = the particles positions
% _________________________________________________________________________
% OUTPUT ARGUMENTS
%
% v         = random directions
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

N = size(particles,1);
% Generate a random pair of indices for each sample
I = zeros(2,N);                         % Create storage for N pairs of indexes
I(1,:) = randi(N-1,1,N);                % For each seed, select a random integer between 1 and (N-1)
I(2,:) = randi(N-2,1,N);                % For each seed, select a random integer between 1 and (N-2)
I(2,:) = I(2,:) + (I(2,:) >=I(1,:));    % Adjust the second random integer if it coincides with the first random integer
I = I +(I>=[(1:N);(1:N)]);              % Adjust both integers to prevent either of them being the same with its own sample index.

% Generate random vector v for each sample
v = 2*(particles(I(2,:),:) - particles(I(1,:),:));

end