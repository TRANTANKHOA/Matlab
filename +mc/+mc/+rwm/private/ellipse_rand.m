function z = ellipse_rand(target)
%ellipse_rand 
%
% _________________________________________________________________________
% INPUT ARGUMENTS
%
% _________________________________________________________________________
% OUTPUT ARGUMENTS
%
% -------------------- Copyright (C) 2016 Khoa T. Tran --------------------
%   Created:        08-Aug-2016
%   Last edited:    08-Aug-2016
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

%% Generating uniform random vector on a sphere
z = randn(target.M,target.Dims);
z = bsxfun(@rdivide,z,sqrt(sum(z.^2,2)));
z = bsxfun(@times,z,rand(target.M,1).^(1/target.Dims));
z = z*target.chol*target.s;