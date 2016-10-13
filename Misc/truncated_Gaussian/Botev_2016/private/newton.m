function x=newton(p,l,u)
q=@(x)(erfcx(x/sqrt(2))/2); % define Q function
ql=q(l);qu=q(u);
l=l.^2;u=u.^2;
% set initial value for Newton iteration
x=sqrt(l-2*reallog(1+p.*expm1(l/2-u/2)));
% initialize Newton method
err=Inf;
while err>10^-10
    del=-q(x)+(1-p).*exp(.5*(x.^2-l)).*ql+p.*exp(.5*(x.^2-u)).*qu;
    x=x-del; % Newton's step
    err=max(abs(del)); % find the maximum error
end