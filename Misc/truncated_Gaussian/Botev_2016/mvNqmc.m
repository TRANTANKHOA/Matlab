function  est=mvNqmc(l,u,Sig,n)
%% truncated multivariate normal cumulative distribution (qmc version)
% computes an estimator of the probability Pr(l<X<u),
% where 'X' is a zero-mean multivariate normal vector
% with covariance matrix 'Sig', that is, X~N(0,Sig)
% infinite values for vectors 'u' and 'l' are accepted;
% 
% This version uses a Quasi Monte Carlo (QMC) pointset
% of size ceil(n/12) and estimates the relative error
% using 12 independent randomized QMC estimators; QMC
% is slower than ordinary Monte Carlo (see my mvncdf.m), 
% but is also likely to be more accurate when d<50.
%
% output:      structure 'est' with
%              1. estimated value of probability Pr(l<X<u)
%              2. estimated relative error of estimator
%              3. theoretical upper bound on true Pr(l<X<u)
%              
% * Remark: If you want to estimate Pr(l<Y<u),
%           where Y~N(m,Sig) has mean vector 'm',
%           then use 'mvNqmc(Sig,l-m,u-m,n)'.
%
% * Example:
%  clear all,clc,d=25;
%  l=ones(d,1)*5;u=Inf(d,1);
%  Sig=0.5*eye(d)+.5*ones(d,d);
%  est=mvNqmc(l,u,Sig,10^4) % output of our method
% % Executing Matlab's toolbox\stats\stats\mvncdf.m
% % with n=10^7 below is slow and inaccurate
%  options=optimset('TolFun',0,'MaxFunEvals',10^7,'Display','iter');
%  [prob,err]=mvncdf(l,u,zeros(d,1),Sig,options)
%
% See also: mvncdf, mvrandn
%
% For more help, see <a href="matlab:
% doc">Truncated Multivariate Student & Normal</a> documentation at the bottom.

% Reference: Z. I. Botev (2017), _The Normal Law Under Linear Restrictions:
% Simulation and Estimation via Minimax Tilting_, Journal of the Royal 
% Statistical Society, Series B, Volume 79, Part 1, pp. 1-24
l=l(:); u=u(:); % set to column vectors
d=length(l); % basic input check
if  (length(u)~=d)|(d~=sqrt(prod(size(Sig)))|any(l>u))
    error('l, u, and Sig have to match in dimension with u>l')
end
% Cholesky decomposition of matrix
[L, l, u]=cholperm( Sig, l, u ); D=diag(L);
if any(D<eps)
    warning('Method may fail as covariance matrix is singular!')
end
L=L./repmat(D,1,d);u=u./D; l=l./D; % rescale
L=L-eye(d); % remove diagonal
% find optimal tilting parameter via non-linear equation solver
options=optimset('Diagnostics','off','Display','off',...
    'Algorithm','trust-region-dogleg','Jacobian','on');
[soln,fval,exitflag] = fsolve(@(x)gradpsi(x,L,l,u),zeros(2*(d-1),1),options);
if exitflag~=1
    warning('Method may fail as covariance matrix is close to singular!')
end
x=soln(1:(d-1));mu=soln(d:(2*d-2)); % assign saddlepoint x* and mu*
for i=1:12 % repeat randomized QMC
    p(i)=mvnprqmc(ceil(n/12),L,l,u,mu);
end
est.prob=mean(p); % average of QMC estimates
est.relErr=std(p)/sqrt(12)/est.prob; % relative error
est.upbnd=psy(x,L,l,u,mu); % compute psi star
if est.upbnd<-743
    warning('Natural log of probability is less than -743, yielding 0 after exponentiation!')
end
est.upbnd=exp(est.upbnd);


