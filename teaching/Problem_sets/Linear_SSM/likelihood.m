function l = likelihood(sys,theta)
%likelihood A script to compute the likelihood
%
% _________________________________________________________________________
% INPUT ARGUMENTS
% sys   = state space model
% theta = parameter of interest
% _________________________________________________________________________
% OUTPUT ARGUMENTS
% l     = log-likelihood of theta
% -------------------- Copyright (C) 2016 Khoa T. Tran --------------------
%   Created:        09-Aug-2016
%   Last edited:    09-Aug-2016
%   Matlab version: 8.6.0.267246 (R2015b)
%   Email:          khoa.tran@unsw.edu.au
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
sys.sigma   = sqrt(exp(theta(1)));
sys.tau     = sqrt(exp(theta(2)));
sys.k       = sys.tau/sqrt(1-sys.phi^2);

sys.Va      = sys.tau^2*eye(sys.T); sys.Va(1,1) = 0;
sys.Ainv_g  = sys.A\sys.g;
sys.Vx      = sys.Ainv_g*sys.Ainv_g'*sys.k^2 + sys.Va;
sys.Vy      = sys.Vx + sys.sigma^2*eye(sys.T);
[L,R]       = ldl(sys.Vy); R = diag(R);
eps         = L\sys.y;
l           = 0.5*(sum(eps.^2./R) + sum(log(R)));
