function particles = my_mvnrnd(target)
%MY_MVNRND Creating M multivariate Gaussian random row vectors
%
% _________________________________________________________________________
% INPUT ARGUMENTS
%
% target.mu     = mean
% target.chol   = Cholesky factor of the covariance matrix
% target.M      = desired number of samples
% target.Dims   = dimensionality
% _________________________________________________________________________
% OUTPUT ARGUMENTS
% 
% particles     = multivariate Gaussian random vectors
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

particles = bsxfun(@plus,randn(target.M,target.Dims) * target.chol,target.mu);

end
    