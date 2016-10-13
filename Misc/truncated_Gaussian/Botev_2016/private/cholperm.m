function [ L, l, u, perm ] = cholperm( Sig, l, u )
%  Computes permuted lower Cholesky factor L for Sig
%  by permuting integration limit vectors l and u.
%  Outputs perm, such that Sig(perm,perm)=L*L'.
%
% Reference: 
%  Gibson G. J., Glasbey C. A., Elston D. A. (1994),
%  "Monte Carlo evaluation of multivariate normal integrals and
%  sensitivity to variate ordering", 
%  In: Advances in Numerical Methods and Applications, pages 120--126

d=length(l);perm=1:d; % keep track of permutation
L=zeros(d,d);z=zeros(d,1);
for j=1:d
    pr=Inf(d,1); % compute marginal prob.
    I=j:d; % search remaining dimensions
    D=diag(Sig);
    s=D(I)-sum(L(I,1:j-1).^2,2);
    s(s<0)=eps;s=sqrt(s);
    tl=(l(I)-L(I,1:j-1)*z(1:j-1))./s;
    tu=(u(I)-L(I,1:j-1)*z(1:j-1))./s;
    pr(I)=lnNpr(tl,tu);
    % find smallest marginal dimension
    [dummy,k]=min(pr);
    % flip dimensions k-->j
    jk=[j,k];kj=[k,j];
    Sig(jk,:)=Sig(kj,:);Sig(:,jk)=Sig(:,kj); % update rows and cols of Sig
    L(jk,:)=L(kj,:); % update only rows of L
    l(jk)=l(kj);u(jk)=u(kj); % update integration limits
    perm(jk)=perm(kj); % keep track of permutation
    % construct L sequentially via Cholesky computation
    s=Sig(j,j)-sum(L(j,1:j-1).^2);
    if s<-0.01
        error('Sigma is not positive semi-definite')
    end
    s(s<0)=eps;L(j,j)=sqrt(s);
    L(j+1:d,j)=(Sig(j+1:d,j)-L(j+1:d,1:j-1)*(L(j,1:j-1))')/L(j,j);
    % find mean value, z(j), of truncated normal:
    tl=(l(j)-L(j,1:j-1)*z(1:j-1))/L(j,j);
    tu=(u(j)-L(j,1:j-1)*z(1:j-1))/L(j,j);
    w=lnNpr(tl,tu); % aids in computing expected value of trunc. normal
    z(j)=(exp(-.5*tl.^2-w)-exp(-.5*tu.^2-w))/sqrt(2*pi);
end