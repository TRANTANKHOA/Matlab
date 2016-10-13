function y = log_normpdf(x,SIGMA,dim)
%NORMPDF Natural logarirhm of the normal probability density function (pdf).
%   Y = LOG_NORMPDF(X,SIGMA) returns the log pdf of the normal distribution
%   with mean 0 and variance SIGMA, evaluated at the values in X.
%
%   X and SIGMA are matched by 'bsxfun' and by default every row of x is
%   considered as one observation from a spherical Gaussian distribution
%   with diagonal covariance matrix. See code for more specifics.
%
%   See also NORMCDF, NORMFIT, NORMINV, NORMLIKE, NORMRND, NORMSTAT.

%   References:
%      [1] Evans, M., Hastings, N., and Peacock, B. (1993) Statistical
%          Distributions, 2nd ed., Wiley, 170pp.

%   Copyright 1993-2004 The MathWorks, Inc. Modified by Khoa Tran, 2015.


if nargin<1
    error(message('stats:normpdf:TooFewInputs'));
elseif nargin < 2
    SIGMA = 1;dim = 2;
elseif nargin <3
    dim = 2;
end

% Return NaN for out of range parameters.
SIGMA(SIGMA <= 0) = NaN;

% Below is basically a sum of all the independent logdensities over all
% dimensions. There is no calculation for det(SIGMA), so this fucntion is
% not suitable for the general multivariate Gaussian density.
y = -0.5 * sum( bsxfun( @plus , bsxfun( @rdivide, x.^2 , SIGMA) , log(SIGMA) + log(2*pi)) ,dim);

