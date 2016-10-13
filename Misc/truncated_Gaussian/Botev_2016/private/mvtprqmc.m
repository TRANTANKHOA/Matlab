function p=mvtprqmc(n,L,l,u,nu,mu)
% computes P(l<X<u), where X is student with
% Sig=L*L', zero mean vector, and degrees of freedom 'nu';
% exponential tilting uses parameter 'mu';
% Quasi Monte Carlo uses 'n' samples;
d=length(l); % Initialization
eta=mu(d); 
Z=zeros(d,n); % create array for variables
% QMC pointset
x = scramble(sobolset(d),'MatousekAffineOwen');
% precompute constants 
const=log(2*pi)/2-gammaln(nu/2)-(nu/2-1)*log(2)+lnNpr(-eta,Inf)+.5*eta^2;
R=eta+norminvp(x(1:n,d),-eta(ones(n,1)),Inf(n,1))'; % simulate R~N(eta,1) with R>0
p=(nu-1)*log(R)-eta*R+const; % compute Likelihood Ratio for R
R=R/sqrt(nu); % scale parameter divided by nu
for k=1:(d-1)
    % compute matrix multiplication L*Z
    col=L(k,1:k)*Z(1:k,:);
    % compute limits of truncation
    tl=R*l(k)-mu(k)-col;
    tu=R*u(k)-mu(k)-col;
    %simulate N(mu,1) conditional on [tl,tu] via QMC
    Z(k,:)=mu(k)+norminvp(x(1:n,k),tl,tu);
    % update likelihood ratio
    p = p+lnNpr(tl,tu)+.5*mu(k)^2-mu(k)*Z(k,:);
end
% deal with final Z(d) which need not be simulated
col=L(d,:)*Z;tl=R*l(d)-col;tu=R*u(d)-col;
p=p+lnNpr(tl,tu); % update LR corresponding to Z(d)
p=mean(exp(p)); % now switch back from logarithmic scale