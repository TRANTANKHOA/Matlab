function x=tn(l,u)
% samples a column vector of length=length(l)=length(u)
% from the standard multivariate normal distribution,
% truncated over the region [l,u], where -a<l<u<a for some
% 'a' and l and u are column vectors;
% uses acceptance rejection and inverse-transform method;
tol=2; % controls switch between methods
% threshold can be tuned for maximum speed for each platform
% case: abs(u-l)>tol, uses accept-reject from randn
I=abs(u-l)>tol; x=l;
if any(I)
    tl=l(I); tu=u(I); x(I)=trnd(tl,tu);
end
% case: abs(u-l)<tol, uses inverse-transform
I=~I;
if any(I)
    tl=l(I); tu=u(I); pl=erfc(tl/sqrt(2))/2; pu=erfc(tu/sqrt(2))/2;
    x(I)=sqrt(2)*erfcinv(2*(pl-(pl-pu).*rand(size(tl))));
end
