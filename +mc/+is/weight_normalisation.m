function [w, ESS] = weight_normalisation( log_w )
%WEIGHT_NORMALISATION Converting the weights from the log-scale to
% linear-scale and normalising them to sum to one.
% Calculating the effective sample size and give an error when the sum of 
% weights is zero.
%
% _________________________________________________________________________
% INPUT ARGUMENTS
%
% log_w = logarithm of the weights
% _________________________________________________________________________
% OUTPUT ARGUMENTS
%
% w     = the weights
% ESS   = effective sample sizes
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


log_sum = lse(log_w);               % Summing all the weights on the log scale

if log_sum ~= -Inf                  % Checking of the sum of weights is not 0.
    log_w = log_w-log_sum;          % Normalising the weights in the log scale
    w = exp(log_w);                 % Convert to linear scale
    w = w/sum(w);                   % Normalising again
    ESS = exp(-lse(log_w*2));       % ESS = 1/sum(w.^2)
else
    error('Sum of weights is zero');   % Zero sum of weights means the importance density is completely incorrect.
end

end

