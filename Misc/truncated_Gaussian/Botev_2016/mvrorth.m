function  rv=mvrorth(l,u,Sig,n)
%% multivariate normal generator suitable for orthant regions
% simulates 'n' random vectors exactly/perfectly distributed
% from the d-dimensional N(0,Sig) distribution (zero-mean normal
% with covariance 'Sig') conditional on l<X<u;
% method suitable for orthant regions such as when 'l=0' and 'u=Inf';
% output:   'd' times 'n' array 'rv' storing random vectors;
%
% * Example: See toolbox documentation for Probit Bayesian Application.
% * Notes: Algorithm may not work if 'Sig' is close to being rank deficient.
%
% See also: mvNcdf, mvNqmc, mvrandn
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
% Cholesky decomposition of matrix with permuation
[Lfull,l,u,perm]=cholorth(Sig,l,u); % outputs the permutation
D=diag(Lfull);L=Lfull./repmat(D,1,d);u=u./D; l=l./D; % rescale
L=L-eye(d); % remove diagonal
% find optimal tilting parameter non-linear equation solver
options=optimset('Diagnostics','off','Display','off',...
    'Algorithm','trust-region-dogleg','Jacobian','on');
soln= fsolve(@(x)gradpsi(x,L,l,u),zeros(2*(d-1),1),options);
x=soln(1:(d-1));mu=soln(d:(2*d-2));
% compute psi star
psistar=psy(x,L,l,u,mu);
% start acceptance rejection sampling
rv=[]; accept=0; iter=0;
while accept<n % while # of accepted is less than n
    [logpr,Z]=mvnrnd(n,L,l,u,mu); % simulate n proposals
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
rv=rv(:,1:n); % cut-down the array to desired n samples
rv=Lfull*rv; % reverse scaling of L
rv=rv(order,:); % reverse the Cholesky permutation



