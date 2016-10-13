function [ Lp, lp, up, perm] = cholorth( Sig, l, u )
%  Computes permuted lower Cholesky factor L for Sig
%  by permuting integration limit vectors l and u, suitable for
%  the computation of orthant probabilities as in the Probit Bayesian Regression;
%  Outputs 'perm', such that Sig(perm,perm)=Lp*Lp'.
d=length(l);Lfull=chol(Sig,'lower');D=diag(Lfull);
if any(D<eps)
    warning('Method may fail as covariance matrix is singular!')
end
L=Lfull./repmat(D,1,d);up=u./D; lp=l./D; % rescale
L=L-eye(d); % remove diagonal
% find optimal tilting parameter non-linear equation solver
options=optimset('Diagnostics','off','Display','off',...
    'Algorithm','trust-region-dogleg','Jacobian','on');
[soln,fval,exitflag] = fsolve(@(x)gradpsi(x,L,l,u),zeros(2*(d-1),1),options);
if exitflag~=1
    warning('Method may fail as covariance matrix is close to singular!')
end
x=soln(1:(d-1));x(d)=0;x=x(:); % adjust for case of 2-dimensions
[dum,perm]=sort(Lfull*x,'ascend'); % find permutation
Lp=chol(Sig(perm,perm),'lower'); % compute new Cholesky factor
lp=l(perm); up=u(perm); % permute limits
