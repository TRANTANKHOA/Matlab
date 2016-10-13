function  [Z,R]=tregress(l,u,Sig,df,n)
%% truncated student generator for Bayesian regression simulation
% simulates 'n' random vectors exactly/perfectly distributed
% from the d-dimensional t_nu(0,Sig) distribution (zero-mean student
% with scale matrix 'Sig' and deg. of freedom df), conditional on l<X<u,
% by outputting a vector 'Z' and scale 'R' such that sqrt(df)*Z/R has the
% desired truncated student pdf; Inf value for 'l' and 'u' accepted;
% output:   1. 'd' times 'n' array 'Z' storing random vectors
%           2. '1' times 'n' scale vector 'R' so that sqrt(df)*Z/R follows
%                                                the truncated student pdf
%
% * Notes: Algorithm may not work if 'Sig' is close to being rank deficient.
%
% See also: mvrandt, mvNcdf, mvNqmc, mvTcdf, mvTqmc, mvrorth
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
% Cholesky decomposition of matrix with permuation
[Lfull,l,u,perm]=cholperm(Sig,l,u); % outputs the permutation
D=diag(Lfull);
if any(D<eps)
    warning('Method may fail as covariance matrix is singular!')
end
L=Lfull./repmat(D,1,d);u=u./D; l=l./D; % rescale
L=L-eye(d); % remove diagonal
% find optimal tilting parameter non-linear equation solver
options=optimset('Diagnostics','off','Display','off',...
    'Algorithm','trust-region-dogleg');
xo=zeros(2*d,1); xo(2*d)=sqrt(df); xo(d)=log(sqrt(df));
[soln,fval,exitflag] = fsolve(@(x)gradpsiT(x,L,l,u,df),xo,options);
if exitflag~=1
    warning('Method may fail as covariance matrix is close to singular!')
end
% assign saddlepoint x* and mu*
soln(d)=exp(soln(d));x=soln(1:d);mu=soln((d+1):end);
% compute psi star
psistar=psyT(x,L,l,u,df,mu);
% start acceptance rejection sampling
rv=[]; accept=0; iter=0;
while accept<n % while # of accepted is less than n
    [logpr,Z,R]=mvtrnd(n,L,l,u,df,mu); % simulate n proposals
    Z=[Z;R]; % store the scale temporarily in array Z
    idx=-log(rand(1,n))>(psistar-logpr); % acceptance tests
    rv=[rv,Z(:,idx)];  % accumulate accepted
    accept=size(rv,2); % keep track of # of accepted
    iter=iter+1;  % keep track of while loop iterations
    if iter==10^3 % if iterations are getting large, give warning
        warning('Acceptance prob. smaller than 0.001')
    elseif iter>10^4 % if iterations too large, seek approximation only
        accept=n;rv=[rv,Z]; % add the approximate samples
        warning('Sample is only approximately distributed.')
    end
end
% finish sampling; postprocessing
[dum,order]=sort(perm,'ascend');
Z=Z(:,1:n); % cut-down the array to desired n samples
R=Z(d+1,:); Z=Z(1:d,:); % unpack the scale from 'Z'
Z=Lfull*Z; % reverse scaling of L
Z=Z(order,:); % reverse the Cholesky permutation




