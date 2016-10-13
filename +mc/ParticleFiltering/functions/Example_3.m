%% LINEAR GAUSSIAN STATE SPACE MODEL
% Given the model with 1 state and 1 parameter
% 
% $$x_{t+1} = 0.7x_t+v_t, v_t\sim{Normal}(0,10^\theta)$$
%
% $$y_t = 0.5x_t+e_t, e_t\sim{Normal}(0,10^{-1})$$
%
% $$\theta\sim{Uniform}(R)$$
% 
% The true parameter is $\theta_0 = 0$, which is used in the next section
% to generate the record $Y_N = y_{1:N}$ for $N=100$.

%% Simulating the output record
clear;close all;clc
N = 1e2;        % Generate 100 data points for simulation
theta_0 = 0;    % Define the true parameter
%%
% The initial state is generated from its stationary distribution, which is
% derived to be as follows
sigma_x = 1/(1-0.7^2);sigma_y = (0.1);
x = randn*sqrt(sigma_x)*ones(1,1,N);
y = x;          % Allocating memory for the output
for t=2:N       % Simulating the state trajectory
    x(t) = x(t-1)*0.7+randn;
end
noise = sqrt(sigma_y)*randn(size(y));
y = x*0.5+noise;% Simulating the output trajectory
% Printing the state vs output trajectories
f1 = figure;
[hAx,hLine1,hLine2] = plotyy(1:N,x(:),1:N,y(:),'plot','stem');
title('Simulated Output and State Trajectory')
xlabel('Time index')
ylabel(hAx(1),'State') % left y-axis
ylabel(hAx(2),'Output') % right y-axis
hLine1.LineStyle = '-.';
hLine1.Marker = '*';
% hLine2.LineStyle = ':';
set(hAx(2),'Ylim',get(hAx(1),'Ylim')/2)
%% Solving p(X|theta,Y)
% Draw m values of theta, then draw m values of x for each theta
m = 10;log_m = log(m);
theta = randn(m,1);     
X = repmat(bsxfun(@times,randn(m,m,1),sqrt(10.^theta./(1-0.7^2))),1,1,N);
% m(theta) x N(steps) x m(x)
for t=2:N       % Simulating the state trajectory
    X(:,:,t) = X(:,:,t-1)*0.7+bsxfun(@times,randn(m,m,1),sqrt(10.^theta));
end%% 
% Printing the state trajectory to compare with the conditional mean
% estimate later
f2 = figure; hold on
plot(x(:))
qt = quantile(X,[5 95]/100,2);
ln1 = line('XData',(1:t)','YData',mean(qt(:,1,:),1),'LineStyle','-.','Color','g');
ln2 = line('XData',(1:t)','YData',mean(qt(:,2,:),1),'LineStyle','-.','Color','m');
xlabel('Time index')
ylabel('State') 
legend('True states','5^{th} percentile','95^{th} percentile')
drawnow
%% Particle Filtering Operations 
% We can sequentially migrate the parameter-state particles through the
% filtering distribution sequence $p(\theta,x_t | y_{1:t})$. Ocasionally,
% the algorithm will automatically create intermediate densities between
% posterior and prior densities if the "gap" is too large. We demand that
% each consecutive density cannot be too different from its predecessor
% such that at least 75% of the particles remain "important"
% 
log_rho = log(0.95); ESS_0 = log_rho + log_m;
%
% The particles will then be propagated through this sequence of densities
% by importance-sampling-resampling and then slice sampling. We perform
% calculation on the log-scale throughout this example.
% 
% Allocating memory for the loops
X1 = X;   X2 = X(:,:,1);    X3 = X2;    X4 = X3;    e = X;
    
%% Define the prior and posterior

prior = log_prior_theta(X,theta);
posterior = log_likelihood_theta(X,y)+prior;

% The first transition density is coincided with the prior.
transit = prior;            
lambda_0 = zeros(m,1);lambda_1 = lambda_0;  % The transition index is set to be zero 
flags = lambda_0<1;
while any(flags)            % The loop stops when $\lambda = 0$
    %% Define next transition index $\lambda_1$
    % Redefining the transition density if the previous search for
    % $\lambda_1$ was failed.
    chk = lambda_0(flags) ~= lambda_1(flags);
    if any(chk)     
        transit(chk,:,1) = lwse_3d([posterior(chk,:,1) prior(chk,:,1)],[lambda_0(chk) (1-lambda_0(chk))]);
    end
    lambda_1(flags) = 1; % Start with the most ambitious increment
    % Performing a bisection search for $\lambda$
    bis_flg = flags;
    while any(bis_flg)
        delta(bis_flg,1) = lambda_1(bis_flg)-lambda_0(bis_flg);
        P1(bis_flg,:,1) = lwse_3d(cat(3, posterior(bis_flg,:,1), prior(bis_flg,:,1)),cat(3, lambda_1(bis_flg), 1-lambda_1(bis_flg)),3);
        log_omega(bis_flg,:,1) = P1(bis_flg,:,1)-transit(bis_flg,:,1);
        % Flagging particles with infinite weight, if any. Then replacing them
        % with the next highest weight in the pool.
        Inf_flg(bis_flg,:,1) = log_omega(bis_flg,:,1)==Inf;
        if any(Inf_flg(:))
            log_omega(Inf_flg) = max(log_omega(~Inf_flg))*ones(size(log_omega(Inf_flg)));
        end
        % Calculating the effective sample size
        
        log_sum(bis_flg,1) = lse_3d(log_omega(bis_flg,:,1),2);
        sum_flg(bis_flg,1) = log_sum(bis_flg,1) ~= -Inf ;   % Checking of the sum of weights is not 0.
        log_omega(sum_flg,:,1) = bsxfun(@minus,log_omega(sum_flg,:,1),log_sum(sum_flg));   % Normalising the weights
        log_ESS(sum_flg,1) = -lse_3d(log_omega(sum_flg,:,1)*2,2);            % ESS = 1/sum(w.^2)
        log_ESS(~sum_flg,1) = -Inf;         % Zero effective sample size
        
        
        % Evaluate the condition for breaking out of the loop
        ESS_flg(bis_flg,1) = log_ESS(bis_flg,1)  >= ESS_0;
        lambda_0(ESS_flg & bis_flg) = lambda_1(ESS_flg & bis_flg);
        % Sometime even a very small increment in $\lambda$ is already
        % too big for the predefined threshold. In such case, we keep
        % lambda_0 the same and use the slice sampler to predistribute
        % the pool of particles toward lambda_1. This is also why we
        % had to guard our code against particles with infinite weight.
        bis_flg = ~(ESS_flg | (delta<1e-6));
        lambda_1(bis_flg) = (lambda_1(bis_flg)+lambda_0(bis_flg))/2; 
    end
    clc
    formatSpec = 'Time index is t = %d.\n';
    fprintf(formatSpec,t)
    formatSpec = 'Effective sample size is %2.1f %% samples.\n';
    fprintf(formatSpec,median(exp(log_ESS))/m*100)
    formatSpec = 'Lambda = %1.6f and delta = %1.6f.\n';
    fprintf(formatSpec,median([lambda_1 delta]))
    %% Importance Resampling
    % Multinomial sampling is more accurate, yet slower
    % 
    MUL_P = exp(log_omega); MUL_P = bsxfun(@rdivide,MUL_P,sum(MUL_P,2));
    r = mnrnd(m,MUL_P); mark_end = zeros(m,1); mark_begin = mark_end;
    rspl = zeros(m,m);
    for i=1:m
        r_flg = r(:,i)~=0;
        mark_begin(r_flg,1) = mark_end(r_flg,1)+1;
        mark_end(r_flg,1) = mark_begin(r_flg,1)+r(r_flg,i)-1;
        rw = nonzeros((1:m).*r_flg');
        for j = 1:length(rw)
            rspl(rw(j),mark_begin(rw(j)):mark_end(rw(j))) = i;
        end
    end
    for i=1:m
        for j=1:m
            X1(i,j,:) = X(i,rspl(i,j),:);
            X2(i,j,:) = prior(i,rspl(i,j),:);
            X3(i,j,:) = posterior(i,rspl(i,j),:);
            X4(i,j,:) = P1(i,rspl(i,j),:);
        end
    end      
    X = X1; prior = X2; posterior = X3; P1 = X4;
    % 
    % Instead, we can also simply throw away the lower 25th percentile and
    % replace them with the upper 25th percentile as follows
    %
    %     [~,Idx] = sort(log_omega);
    %     X(Idx(1:round(m/4)),:) = X(Idx(m-round(m/4)+1:m),:);
    %     prior(Idx(1:round(m/4)),:) = prior(Idx(m-round(m/4)+1:m),:);
    %     posterior(Idx(1:round(m/4)),:) = posterior(Idx(m-round(m/4)+1:m),:);
    %     P1(Idx(1:round(m/4)),:) = P1(Idx(m-round(m/4)+1:m),:);
    %
%% Adative Direction Slice sampling - One step only
    m_w = 100;                  % Maximum number of expansions
    height = P1-exprnd(1,m,m);  % Draw a random height for each seed.
    %% Generate random indices
    I = zeros(m,m,2);
    I(:,:,1) = randi(m-1,m,m);        % For each seed, select a random integer between 1 and (m-1)
    I(:,:,2) = randi(m-2,m,m);        % For each seed, select a random integer between 1 and (m-2)
    I(:,:,2) = I(:,:,2) + (I(:,:,2) >=I(:,:,1));    % Adjust the second random integer if it coincides with the first random integer
    I = I +(I>=repmat(1:m,[m,1,2]));      % Adjust both integers to prevent selecting the 'current' seed
    %% Generate random directions
    for i=1:m
        for j=1:m
            e(i,j,:) = 2.*(X(i,I(i,j,2),:) - X(i,I(i,j,1),:));
            while ~any(e(i,j,:))
                e(i,j,:) = 2.*(X(i,I(i,j,2),:) - X(i,I(i,randi(m),1),:));
            end
        end
    end
    %% Random initial bounds
    b = rand(m,m);
    a = b - 1;
    %% Expansion upper bound
    m_b = randi(m_w+1,m,m)-1;
    for i=1:m
        cnt=0;    slc_flg = true(1,m);
        while any(slc_flg)
            X1(i,slc_flg,:) = X(i,slc_flg,:)+ bsxfun(@times,e(i,slc_flg,:),b(i,slc_flg));
            A1 = log_prior_theta(X1(i,slc_flg,:),theta(i));
            B1 = log_likelihood_theta(X1(i,slc_flg,:),y)+A1;
            slc_flg(slc_flg) = lwse_3d(cat(3,B1,A1),cat(3,lambda_1(i),(1-lambda_1(i))),3)>=height(i,slc_flg);
            slc_flg(slc_flg) = slc_flg(slc_flg) & (cnt<m_b(i,slc_flg));
            b(i,slc_flg) = b(i,slc_flg)+slc_flg(slc_flg);
            cnt = cnt+1;
        end
        formatSpec = 'Expansions on upper bound is %d.\n';
        fprintf(formatSpec,cnt-1)
    end
    %% Expansion lower bound
    m_a = m_w - m_b;
    for i=1:m
        cnt=0;    slc_flg = true(1,m);
        while any(slc_flg)
            X1(i,slc_flg,:) = X(i,slc_flg,:)+ bsxfun(@times,e(i,slc_flg,:),a(i,slc_flg));
            A1 = log_prior_theta(X1(i,slc_flg,:),theta(i));
            B1 = log_likelihood_theta(X1(i,slc_flg,:),y)+A1;
            slc_flg(slc_flg) = lwse_3d(cat(3,B1,A1),cat(3,lambda_1(i),(1-lambda_1(i))),3)>=height(i,slc_flg);
            slc_flg(slc_flg) = slc_flg(slc_flg) & (cnt<m_b(i,slc_flg));
            a(i,slc_flg) = a(i,slc_flg)-slc_flg(slc_flg);
            cnt = cnt+1;
        end
        formatSpec = 'Expansions on lower bound is %d.\n';
        fprintf(formatSpec,cnt-1)
    end
    %% Proposing and shrinking bounds
    for i = 1:m
        slc_flg = true(1,m);
        innovations = zeros(1,m);
        while any(slc_flg)
            innovations(slc_flg) = (b(i,slc_flg)-a(i,slc_flg)).*rand(1,nnz(slc_flg))+a(i,slc_flg);
            X1(i,slc_flg,:) = X(i,slc_flg,:)+bsxfun(@times,e(i,slc_flg,:),innovations(slc_flg)); 
            prior(i,slc_flg) = log_prior_theta(X1(i,slc_flg,:),theta(i));
            posterior(i,slc_flg) = log_likelihood_theta(X1(i,slc_flg,:),y) + prior(i,slc_flg);
            transit(i,slc_flg) = lwse_3d(cat(3,posterior(i,slc_flg), prior(i,slc_flg)),cat(3,lambda_1(i), 1-lambda_1(i)),3);
            slc_flg(slc_flg) = transit(i,slc_flg)<height(i,slc_flg);
            idx = slc_flg;
            idx(slc_flg)  = innovations(slc_flg)>0;
            b(i,idx) = innovations(idx);
            a(i,~idx & slc_flg) = innovations(~idx & slc_flg);
        end
    end
    X  = X1;
figure(f2); hold on
qt = quantile(X,[5 95]/100,2);
ln1.YData = mean(qt(:,1,:),1);
ln2.YData = mean(qt(:,2,:),1);
% xlabel('Time index')
% ylabel('State') 
% legend('True states','5^{th} percentile','95^{th} percentile')
drawnow
end


%% Plotting the precision
% xi = 10.^(-X(:,1));
% figure
% hist(xi,50)