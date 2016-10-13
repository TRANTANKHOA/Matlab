function grad=gradpsiT(y,L,l,u,nu)
% implements gradient of psi(x) to find optimal exponential twisting;
% assumes scaled 'L' with zero diagonal;
d=length(u);c=zeros(d,1);x=c;mu=c;
x(1:(d-1))=y(1:(d-1));mu(1:(d-1))=y(d+1:2*d-1);
r=exp(y(d)); eta=y(2*d); % split along 'mu', 'x', 'r', 'eta'
l=l/sqrt(nu); u=u/sqrt(nu);
% compute now ~l and ~u
c(2:d)=L(2:d,:)*x;lt=r*l-mu-c;ut=r*u-mu-c;
% compute gradients avoiding catastrophic cancellation
w=lnNpr(lt,ut);
pl=exp(-0.5*lt.^2-w)/sqrt(2*pi);
pu=exp(-0.5*ut.^2-w)/sqrt(2*pi);
P=pl-pu;
% output the gradient
dfdx=-mu(1:(d-1))+(P'*L(:,1:(d-1)))';
dfdm= mu-x+P;
l(isinf(l))=0; u(isinf(u))=0;
dfdr=(nu-1)/r-eta+sum(u.*pu-l.*pl);
dfde=eta-r+exp(-0.5*eta^2-lnNpr(-Inf,eta))/sqrt(2*pi);
grad=[dfdx;dfdm(1:d-1);dfdr;dfde];
