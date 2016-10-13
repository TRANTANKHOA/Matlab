function p=mvnprqmc(n,L,l,u,mu)
% computes P(l<X<u), where X is normal with
% 'Cov(X)=L*L' and zero mean vector;
% exponential tilting uses parameter 'mu';
% Quasi Monte Carlo uses 'n' samples;
d=length(l); % Initialization
mu(d)=0;
Z=zeros(d,n); % create array for variables
% QMC pointset
x = scramble(sobolset(d-1),'MatousekAffineOwen');
p=0;
for k=1:(d-1)
    % compute matrix multiplication L*Z
    col=L(k,1:k)*Z(1:k,:);
    % compute limits of truncation
    tl=l(k)-mu(k)-col;
    tu=u(k)-mu(k)-col;
    %simulate N(mu,1) conditional on [tl,tu] via QMC
    Z(k,:)=mu(k)+norminvp(x(1:n,k),tl,tu);
    % update likelihood ratio
    p = p+lnNpr(tl,tu)+.5*mu(k)^2-mu(k)*Z(k,:);
end
% deal with final Z(d) which need not be simulated
col=L(d,:)*Z;tl=l(d)-col;tu=u(d)-col;
p=p+lnNpr(tl,tu); % update LR corresponding to Z(d)
p=mean(exp(p)); % now switch back from logarithmic scale