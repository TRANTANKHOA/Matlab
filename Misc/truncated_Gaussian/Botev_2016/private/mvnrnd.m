function [p,Z]=mvnrnd(n,L,l,u,mu)
% generates the proposals from the exponentially tilted 
% sequential importance sampling pdf;
% output:    'p', log-likelihood of sample
%             Z, random sample 
d=length(l); % Initialization
mu(d)=0;
Z=zeros(d,n); % create array for variables
p=0;
for k=1:d
    % compute matrix multiplication L*Z
    col=L(k,1:k)*Z(1:k,:);
    % compute limits of truncation
    tl=l(k)-mu(k)-col;
    tu=u(k)-mu(k)-col;
    %simulate N(mu,1) conditional on [tl,tu]
    Z(k,:)=mu(k)+trandn(tl',tu');
    % update likelihood ratio
    p = p+lnNpr(tl,tu)+.5*mu(k)^2-mu(k)*Z(k,:);
end