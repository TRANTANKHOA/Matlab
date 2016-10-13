function h = plot3D(x,y,func)
%plot3D Plotting a 3D graph of a 2D function
%
% _________________________________________________________________________
% INPUT ARGUMENTS
% func = function handle
% x,y  = coordinates
% _________________________________________________________________________
% OUTPUT ARGUMENTS
% h    = figure handle
% -------------------- Copyright (C) 2016 Khoa T. Tran --------------------
%   Created:        09-Aug-2016
%   Last edited:    09-Aug-2016
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
x = x(:)';y = y(:)';    % reformat the coordinate
[x,y] = meshgrid(x,y);  % all combinations of x, y
n = length(x)^2;        % number of function evaluations
parfor i = 1:n          % requires parallel computing toolbox
    z(i) = func([x(i) y(i)]);
end
z = reshape(z,size(x)); % put into same size as X, Y
h = meshc(x,y,z);       % contour plot
