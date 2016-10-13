function [p,Z,R]=mvtrnd(n,L,l,u,nu,mu)
% generates the proposals from the exponentially tilted 
% sequential importance sampling pdf;
% output:    'p', log-likelihood of sample
%             Z, Gaussian sample 
%             R, random scale parameter so that sqrt(nu)*Z/R is student

d=length(l); % Initialization
eta=mu(d); mu(d)=0;
Z=zeros(d,n); % create array for variables
% precompute constants
const=log(2*pi)/2-gammaln(nu/2)-(nu/2-1)*log(2)+lnNpr(-eta,Inf)+.5*eta^2;
R=eta+trandn(-eta(ones(n,1)),Inf(n,1))'; % simulate R~N(eta,1) with R>0
p=(nu-1)*log(R)-eta*R+const; % compute Likelihood Ratio for R
%R=R/sqrt(nu); % scale parameter divided by nu
for k=1:d
    % compute matrix multiplication L*Z
    col=L(k,1:k)*Z(1:k,:);
    % compute limits of truncation
    tl=R*l(k)/sqrt(nu)-mu(k)-col;
    tu=R*u(k)/sqrt(nu)-mu(k)-col;
    %simulate N(mu,1) conditional on [tl,tu]
    Z(k,:)=mu(k)+trandn(tl(:),tu(:));
    % update likelihood ratio
    p = p+lnNpr(tl,tu)+.5*mu(k)^2-mu(k)*Z(k,:);
end
