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
N = 1e3;        % Generate 100 data points for simulation
theta_0 = 0;    % Define the true parameter
%%
% The initial state is generated from its stationary distribution, which is
% derived to be as follows
sigma_x = 1/(1-0.7^2);sigma_y = (0.1);
x = randn*sqrt(sigma_x)*ones(N,1);
y = x;          % Allocating memory for the output
for t=2:N       % Simulating the state trajectory
    x(t) = x(t-1)*0.7+randn;
end
noise = sqrt(sigma_y)*randn(size(y));
y = x*0.5+noise;% Simulating the output trajectory
% Printing the state vs output trajectories
f1 = figure;
[hAx,hLine1,hLine2] = plotyy(1:N,x,1:N,y,'plot','stem');
title('Simulated Output and State Trajectory')
xlabel('Time index')
ylabel(hAx(1),'State') % left y-axis
ylabel(hAx(2),'Output') % right y-axis
hLine1.LineStyle = '-.';
hLine1.Marker = '*';
% hLine2.LineStyle = ':';
set(hAx(2),'Ylim',get(hAx(1),'Ylim')/2)
% Printing the state trajectory to compare with the conditional mean
% estimate later
f3 = figure;
plot(x(1:t))
hold on
ln = line('XData',(1:t)','YData',zeros(size(x)));
xlabel('Time index')
ylabel('State') 
legend('True states','Estimated states')
%% Forming the particles at time t = 1
% The initial particles for the state will be generated from its invariant
% distribution as follows
%%
% 
% $$x_1\sim{Normal}(0,\frac{1}{(1-0.7^2)^2})$$
% 
% while the particles for the parameter will be generated from a standard
% Gaussian for convenience.
%%
% 
% $$\theta\sim{Normal}(0,1)$$
% 
m = 1e5;log_m = log(m);
theta = randn(m,1);
state = randn(m,1)/sqrt(1-0.7^2);
%% 
% We form the joint parameter-state particles as follows
X = [theta,state];
Xt = zeros(m,2,N);      % Repository for all X at t = 1:N
% Plotting the bivariate histogram of the current particles
f2 = figure;
[H,C] = hist3(X,[100 100],'LineStyle','none');
H = H/sum(H(:))/(C{1}(2)-C{1}(1))/(C{2}(2)-C{2}(1));
contour(C{1},C{2},H)
title('Initial distribution of $$p(\theta,x_1)$$','Interpreter','Latex')
xlabel('\theta');ylabel('x');
yl = ylim;  xl = xlim;
axis tight

%% Particle Filtering Operations 
% We can sequentially migrate the parameter-state particles through the
% filtering distribution sequence $p(\theta,x_t | y_{1:t})$. Ocasionally,
% the algorithm will automatically create intermediate densities between
% posterior and prior densities if the "gap" is too large. We demand that
% each consecutive density cannot be too different from its predecessor
% such that at least 75% of the particles remain "important"
% 
log_rho = log(0.75);
%
% The particles will then be propagated through this sequence of densities
% by importance-sampling-resampling and then slice sampling. We perform
% calculation on the log-scale throughout this example.
% 
% Allocating memory for the loops
X1 = X;   X2 = X(:,1);    X3 = X2;    X4 = X3;
for t = 1:N
    data = y(t);    scale = 2^(1-t/N);
if t==1
    %%
    % On the first time index, the prior density is available in
    % closed-form. This will no longer be true in subsequent time indexes.
    log_prior = @(X) log_normpdf(X(:,1))+log_normpdf(X(:,2),sigma_x);
    log_posterior = @(X) log_normpdf(bsxfun(@plus,data,-0.5.*X(:,2)),sigma_y)...
        +log_normpdf(X(:,2),sigma_x) + log_normpdf(X(:,1));
    prior = log_prior(X);
    posterior = log_posterior(X);
else
    % From t=2 and onward, the prior density need to be approximate either
    % by a normal density, which is fast and inaccurate, or a kernel
    % density estimator which often has complexity $O(m^2)$ and higher
    % accuracy.
    % 
    % First, we need to simulate the particles for the next state using the
    % state-equation as follows
    X(:,2) = X(:,2).*0.7 + randn(m,1).*sqrt(10.^X(:,1));
    % The first and second moments of the pool of particles are calculated
    mu = mean(X);
    sigma_2 = cov(X);
    % One option is to use the normal density to approximate this prior
    log_prior = @(X) log_mvnpdf(X,mu,sigma_2);
    log_posterior = @(X) log_normpdf(bsxfun(@plus,data,-0.5.*X(:,2)),sigma_y)...
        +log_mvnpdf(X,mu,sigma_2);
    prior = log_prior(X);
    posterior = log_posterior(X);
    % Another option is to use kernel density estimation to be
    % added...here...
end
% The first transition density is coincided with the prior.
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
        if log_ESS - log_m >= log_rho
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
    formatSpec = 'Time index is t = %d.\n';
    fprintf(formatSpec,t)
    formatSpec = 'Effective sample size is %2.1f %% samples.\n';
    fprintf(formatSpec,exp(log_ESS)/m*100)
    formatSpec = 'Lambda = %1.6f and delta = %1.6f.\n';
    fprintf(formatSpec,[lambda_1 delta])
    %% Importance Resampling
    % Multinomial sampling is more accurate, yet slower
    % 
    r = mnrnd(m,exp(log_omega)); mark_end = 0; 
    for i=1:m
        if r(i)~=0
            mark_begin = mark_end+1;mark_end = mark_begin+r(i)-1;
            X1(mark_begin:mark_end,:) = repmat(X(i,:),r(i),1);
            X2(mark_begin:mark_end) = repmat(prior(i),r(i),1);
            X3(mark_begin:mark_end) = repmat(posterior(i),r(i),1);
            X4(mark_begin:mark_end) = repmat(P1(i),r(i),1);
        end
    end 
    X = X1; prior = X2; posterior = X3; P1 = X4;
%     figure(f2)
%     [H,C] = hist3(X1,[100 100],'LineStyle','none'); 
%     H = H/sum(H(:))/(C{1}(2)-C{1}(1))/(C{2}(2)-C{2}(1)); 
%     contour(C{1},C{2},H)
%     title(['Current distribution of $$p(\theta,x_{' num2str(t) '}|Y_{' num2str(t) '})$$'],'Interpreter','Latex')
%     xlabel('\theta');ylabel('x');
%     set(gca,'Xlim',xl*scale,'Ylim',yl*scale) 
%     drawnow
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
    height = P1-exprnd(1,m,1);  % Draw a random height for each seed.
    %% Generate random indices
    I = zeros(2,m);
    I(1,:) = randi(m-1,1,m);        % For each seed, select a random integer between 1 and (m-1)
    I(2,:) = randi(m-2,1,m);        % For each seed, select a random integer between 1 and (m-2)
    I(2,:) = I(2,:) + (I(2,:) >=I(1,:));    % Adjust the second random integer if it coincides with the first random integer
    I = I +(I>=[(1:m);(1:m)]);      % Adjust both integers to prevent selecting the 'current' seed
    %% Generate random directions
    e = 2.*(X(I(2,:),:) - X(I(1,:),:));
    %% Random initial bounds
    b = rand(m,1);
    a = b - 1;
    %% Expansion upper bound
    m_b = randi(m_w+1,m,1)-1;
    i=0;    flags = true(m,1);
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
    i=0;    flags = true(m,1);
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
    flags = true(m,1);
    innovations = zeros(m,1);
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
Xt(:,:,t) = X;              % Recording the particles
ln.YData(t) = mean(X(:,2)); % Updating the plot f3 and then f2
figure(f3)
set(gca,'Xlim',[max(0,t-10) t+10])
drawnow
figure(f2)
[H,C] = hist3(X1,[100 100],'LineStyle','none');
H = H/sum(H(:))/(C{1}(2)-C{1}(1))/(C{2}(2)-C{2}(1));
% mesh(C{1},C{2},H)
contour(C{1},C{2},H)
title(['Current distribution of $$p(\theta,x_{' num2str(t) '}|Y_{' num2str(t) '})$$'],'Interpreter','Latex')
xlabel('\theta');ylabel('x');
set(gca,'Xlim',xl*scale,'Ylim',yl*scale) 
drawnow
end
%% Plotting the precision
xi = 10.^(-X(:,1));
figure
hist(xi,50)