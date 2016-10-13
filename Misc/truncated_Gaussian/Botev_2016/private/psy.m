function p=psy(x,L,l,u,mu)
% implements psi(x,mu); assumes scaled 'L' without diagonal;
d=length(u);x(d)=0;mu(d)=0;x=x(:);mu=mu(:);
% compute now ~l and ~u
c=L*x;l=l-mu-c;u=u-mu-c;
p=sum(lnNpr(l,u)+.5*mu.^2-x.*mu);