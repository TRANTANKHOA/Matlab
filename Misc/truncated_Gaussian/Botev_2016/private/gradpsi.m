function [grad,J]=gradpsi(y,L,l,u)
% implements gradient of psi(x) to find optimal exponential twisting;
% assumes scaled 'L' with zero diagonal;
d=length(u);c=zeros(d,1);x=c;mu=c;
x(1:(d-1))=y(1:(d-1));mu(1:(d-1))=y(d:end);
% compute now ~l and ~u
c(2:d)=L(2:d,:)*x;lt=l-mu-c;ut=u-mu-c;
% compute gradients avoiding catastrophic cancellation
w=lnNpr(lt,ut);
pl=exp(-0.5*lt.^2-w)/sqrt(2*pi);
pu=exp(-0.5*ut.^2-w)/sqrt(2*pi);
P=pl-pu;
% output the gradient
dfdx=-mu(1:(d-1))+(P'*L(:,1:(d-1)))';
dfdm= mu-x+P;
grad=[dfdx;dfdm(1:(d-1))];
if nargout>1 % here compute Jacobian matrix
    lt(isinf(lt))=0; ut(isinf(ut))=0;
    dP=-P.^2+lt.*pl-ut.*pu; % dPdm
    DL=repmat(dP,1,d).*L;
    mx=-eye(d)+DL;
    xx=L'*DL;
    mx=mx(1:d-1,1:d-1);
    xx=xx(1:d-1,1:d-1);
    J=[xx,mx';
        mx,diag(1+dP(1:d-1))];
end