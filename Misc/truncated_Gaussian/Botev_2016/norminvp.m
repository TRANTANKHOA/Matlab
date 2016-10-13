function x=norminvp(p,l,u)
%% normal quantile function with precision
% computes with tail-precision the quantile function
% of the standard normal distribution at 0<=p<=1,
% and truncated to the interval [l,u];
% Inf values for vectors 'l' and 'u' accepted;
%
% * Example 1:
% % Suppose you wish to simulate a random variable
% % 'Z' from the non-standard Gaussian N(m,s^2)
% % conditional on l<Z<u. First compute
%  m=1;l=10;u=20;s=1; 
%  X=norminvp(rand,(l-m)/s,(u-m)/s); % and then set
%  Z=m+s*X
%
% * Example 2:
% % Suppose you desire the median of Z~N(0,1), truncated to Z>9;
%  norminvp(0.5,9,Inf) %our method
% % Matlab's toolbox\stats\stats\norminv.m fails
%   pl=normcdf(9);
%   norminv(0.5*(1-pl)+pl)
%
% See also: trandn
%
% For more help, see <a href="matlab:
% doc">Truncated Multivariate Student & Normal</a> documentation at the bottom.

% Reference: Z. I. Botev (2017), _The Normal Law Under Linear Restrictions:
% Simulation and Estimation via Minimax Tilting_, Journal of the Royal 
% Statistical Society, Series B, Volume 79, Part 1, pp. 1-24

l=l(:);u=u(:);p=p(:); % set to column vectors
if (length(l)~=length(p))|any(l>u)|any(p>1)|any(p<0)
    error('l, u, and p must be the same length with u>l and 0<=p<=1')
end
x=nan(size(l)); % allocate memory
I=(p==1);x(I)=u(I); % extreme values of quantile
J=(p==0);x(J)=l(J);
I=~(I|J); % cases for which 0<x<1
if  any(I)
    x(I)=cases(p(I),l(I),u(I));
end









