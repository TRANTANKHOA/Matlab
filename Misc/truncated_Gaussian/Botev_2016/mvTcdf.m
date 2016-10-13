function  est=mvTcdf(l,u,Sig,df,n)
%% truncated multivariate student cumulative distribution
% computes an estimator of the probability Pr(l<X<u),
% where 'X' is a zero-mean multivariate student vector
% with scale matrix 'Sig' and degrees of freedom 'df', 
% that is, X~t_df(0,Sig) 
% infinite values for vectors 'u' and 'l' are accepted;
% Monte Carlo method uses sample size 'n'; the larger
% the 'n', the smaller the relative error of the estimator;
%      output: structure 'est' with
%              1. estimated value of probability Pr(l<X<u)
%              2. estimated relative error of estimator
%              3. theoretical upper bound on true Pr(l<X<u)
%              Remark: If you want to estimate Pr(l<Y<u),
%                   where Y~t_df(m,Sig) has mean vector 'm',
%                     then use 'mvTcdf(Sig,l-m,u-m,df,n)'.
% * Example:
%  clear all,clc,d=25;
%  l=ones(d,1)/2;u=ones(d,1); df=10;
%  Sig=inv(0.5*eye(d)+.5*ones(d,d));
%  est=mvTcdf(l,u,Sig,df,10^4) % output of our method
% % Executing Matlab's toolbox\stats\stats\mvtcdf.m
% % with n=10^7 below is slow and inaccurate
%  options=optimset('TolFun',0,'MaxFunEvals',10^7,'Display','iter');
%  [prob,err]=mvtcdf(l,u,Sig,df,options)
%
% See also: mvTqmc, mvrandt, mvncdf, mvrandn
%
% For more help, see <a href="matlab:
% doc">Truncated Multivariate Student & Normal</a> documentation at the bottom.

% References:
% [1] Z. I. Botev (2017), _The Normal Law Under Linear Restrictions:
% Simulation and Estimation via Minimax Tilting_, Journal of the Royal 
% Statistical Society, Series B, Volume 79, Part 1, pp. 1-24
%
% [2] Z. I. Botev and P. L'Ecuyer (2015), _EFFICIENT PROBABILITY ESTIMATION 
% AND SIMULATION OF THE TRUNCATED MULTIVARIATE STUDENT-t DISTRIBUTION_, 
% Proceedings of the 2015 Winter Simulation Conference, pages 380-391, 
% (L. Yilmaz, W. Chan, I. Moon, T. Roeder, C. Macal, and M. Rossetti, eds.)

l=l(:); u=u(:); % set to column vectors
d=length(l); % basic input check
if  (d~=sqrt(prod(size(Sig)))|any(l>u))
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
    'Algorithm','trust-region-dogleg');
xo=zeros(2*d,1); xo(2*d)=sqrt(df); xo(d)=log(sqrt(df));
[soln,fval,exitflag] = fsolve(@(x)gradpsiT(x,L,l,u,df),xo,options);
if exitflag~=1
    warning('Method may fail as covariance matrix is close to singular!')
end
% assign saddlepoint x* and mu*
soln(d)=exp(soln(d));x=soln(1:d);mu=soln((d+1):end); 
est=mvtpr(n,L,l,u,df,mu);
% compute psi star
est.upbnd=psyT(x,L,l,u,df,mu); % compute psi star
if est.upbnd<-743
    warning('Natural log of probability is less than -743, yielding 0 after exponentiation!')
end
est.upbnd=exp(est.upbnd);



