function x=normq(p,l,u)
% computes with precision the quantile function
% of the standard normal distribution,
% truncated to the interval [l,u];
% normq assumes 0<l<u and 0<p<1
x=nan(size(l));
I=l>10^5; % if x is too big, no need for Newton iteration
if any(I)
    x(I)=sqrt(l(I).^2-2*reallog(1+p(I).*expm1(l(I).^2/2-u(I).^2/2)));
end
I=~I; % otherwise refine guess via Newton iterations
if any(I)
    x(I)=newton(p(I),l(I),u(I));
end