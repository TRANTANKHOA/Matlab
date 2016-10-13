function [ X ] = sequential_slice_sampling( X, log_prior, log_posterior, rho )
%SEQUENTIAL_SLICE_SAMPLING Given a set of particles with respective prior and
%posterior density function handle, this function will migrate the particle
%from the prior density to the posterior density.
N = size(X,1); log_N = log(N);
d = size(X,2);
X1= zeros(N,d);
log_rho     = log(rho);
prior       = log_prior(X);
posterior   = log_posterior(X);
transit = prior;            
lambda_0 = 0;lambda_1 = 0;  % The transition index is set to be zero 
while lambda_0<1            % The loop stops when $\lambda = 0$
    %% Define next transition index $\lambda_1$
    % Redefining the transition density if the previous search for
    % $\lambda_1$ was failed.
    if lambda_0 ~= lambda_1     
        transit = lwse_row([posterior prior],[lambda_0 (1-lambda_0)]);
    end
    lambda_1 = 1; % Start with the most ambitious increment
    % Performing a bisection search for $\lambda$
    while(true)
        delta = (lambda_1-lambda_0);
        P1 = lwse_row([posterior prior],[lambda_1 (1-lambda_1)]);
        log_omega = P1-transit;
        % Flagging particles with infinite weight, if any. Then replacing them
        % with the next highest weight in the pool.
        Inf_flg = log_omega==Inf;
        if any(Inf_flg)
            log_omega(Inf_flg) = max(log_omega(~Inf_flg))*ones(size(log_omega(Inf_flg)));
        end
        % Calculating the effective sample size
        if lse(log_omega) ~= -Inf   % Checking of the sum of weights is not 0.
            log_omega = log_omega-lse(log_omega);   % Normalising the weights
            log_ESS = -lse(log_omega*2);            % ESS = 1/sum(w.^2)
        else
            log_ESS = -Inf;         % Zero effective sample size
        end
        % Evaluate the condition for breaking out of the loop
        if log_ESS - log_N >= log_rho
            lambda_0 = lambda_1;
            break; 
        elseif delta<1e-6
            % Sometime even a very small increment in $\lambda$ is already
            % too big for the predefined threshold. In such case, we keep
            % lambda_0 the same and use the slice sampler to predistribute
            % the pool of particles toward lambda_1. This is also why we
            % had to guard our code against particles with infinite weight.
            break; 
        else    % Propose the next bisection
            lambda_1 = (lambda_1+lambda_0)/2; 
        end
    end
    clc
    formatSpec = 'Performing Sequential Slice Sampling.\n Effective sample size is %2.1f %% samples.\n';
    fprintf(formatSpec,exp(log_ESS)/N*100)
    formatSpec = 'Lambda = %1.6f and delta = %1.6f.\n';
    fprintf(formatSpec,[lambda_1 delta])
    %% Importance Resampling
    % Performing systematic resampling, which is efficient and often leads
    % to more accurate approximation.
    % 
    w = exp(log_omega);     w = w/sum(w);
    XandP = important_resampling( [X prior posterior P1], w, N );
    X           = XandP(:,1:d);
    prior       = XandP(:,d+1);
    posterior   = XandP(:,d+2);
    P1          = XandP(:,d+3);
    %% Adative Direction Slice sampling - One step only
    m_w = 100;                  % Maximum number of expansions
    height = P1-exprnd(1,N,1);  % Draw a random height for each seed.
    %% Generate random indices
    I = zeros(2,N);
    I(1,:) = randi(N-1,1,N);        % For each seed, select a random integer between 1 and (m-1)
    I(2,:) = randi(N-2,1,N);        % For each seed, select a random integer between 1 and (m-2)
    I(2,:) = I(2,:) + (I(2,:) >=I(1,:));    % Adjust the second random integer if it coincides with the first random integer
    I = I +(I>=[(1:N);(1:N)]);      % Adjust both integers to prevent selecting the 'current' seed
    %% Generate random directions
    e = 2*(X(I(2,:),:) - X(I(1,:),:));
    %% Random initial bounds
    b = rand(N,1);
    a = b - 1;
    %% Expansion upper bound
    m_b = randi(m_w+1,N,1)-1;
    i=0;    flags = true(N,1);
    while any(flags)
        X1(flags,:) = X(flags,:)+bsxfun(@times,e(flags,:),b(flags));
        A1 = log_prior(X1(flags,:));
        B1 = log_posterior(X1(flags,:));
        flags(flags) = lwse_row([B1 A1],[lambda_1 (1-lambda_1)])>=height(flags);
        flags(flags) = flags(flags) & (i<m_b(flags));
        b(flags) = b(flags)+flags(flags);
        i = i+1;
    end
    formatSpec = 'Expansions on upper bound is %d.\n';
    fprintf(formatSpec,i-1)
    %% Expansion lower bound
    m_a = m_w - m_b;
    i=0;    flags = true(N,1);
    while any(flags)
        X1(flags,:) = X(flags,:)+bsxfun(@times,e(flags,:),a(flags));
        A1 = log_prior(X1(flags,:));
        B1 = log_posterior(X1(flags,:));
        flags(flags) = lwse_row([B1 A1],[lambda_1 (1-lambda_1)])>=height(flags);
        flags(flags) = flags(flags) & (i<m_a(flags));
        a(flags) = a(flags)-flags(flags);
        i = i+1;
    end
    formatSpec = 'Expansions on lower bound is %d.\n';
    fprintf(formatSpec,i-1)
    %% Proposing and shrinking bounds
    flags = true(N,1);
    innovations = zeros(N,1);
    while any(flags)
        innovations(flags) = (b(flags)-a(flags)).*rand(nnz(flags),1)+a(flags);
        X1(flags,:) = X(flags,:)+bsxfun(@times,e(flags,:),innovations(flags)); 
        prior(flags) = log_prior(X1(flags,:));
        posterior(flags) = log_posterior(X1(flags,:));
        transit(flags) = lwse_row([posterior(flags) prior(flags)],[lambda_1 (1-lambda_1)]);
        flags(flags) = transit(flags)<height(flags);
        idx = flags;
        idx(flags)  = innovations(flags)>0;
        b(idx) = innovations(idx);
        a(~idx & flags) = innovations(~idx & flags);
    end
    X  = X1;
end

end

