function p=psyT(x,L,l,u,nu,mu)
% implements psi(x,mu); assumes scaled 'L' without diagonal;
d=length(u);r=x(d);eta=mu(d);
x(d)=0;mu(d)=0; x=x(:);mu=mu(:);
l=l/sqrt(nu); u=u/sqrt(nu);
% compute now ~l and ~u
c=L*x;l=r*l-mu-c;u=r*u-mu-c;
p=sum(lnNpr(l,u)+.5*mu.^2-x.*mu);
p=p+log(2*pi)/2-gammaln(nu/2)-(.5*nu-1)*log(2);
p=p+.5*eta^2-r*eta+(nu-1)*reallog(r)+lnNpr(-Inf,eta);