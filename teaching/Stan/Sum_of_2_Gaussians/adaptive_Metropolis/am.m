%   ADAPTIVE METROPOLIS: Generate random samples from a target density phi
%   with a Gaussian proposal. Function syntax is:
%
%   OPT = AM(target,OPT)
%  
%   in which:
%
%   target       Function handle for the target density
%   OPT:
%   OPT.X:       Initial sample
%   OPT.Sig:     Initial covariance matrix to approximate the geometry of 
%                the target distribution
%   OPT.M:       Number of samples to be generated
%   OPT.data     Necessary data to compute the density function  
%   OPT.dsp      Display completion status
%   OPT.hgrd     Number of grid points in histograms, must be multiple of 10
%   OPT.filename Naming the file to store MCMC samples
%   OPT.cov      Adapting the covariance matrix, 1 - adaptive | 0 - fixed
%   OPT.burnin    % Burn-in time
%   OPT.best_X    % Sample with highest joint density
%   OPT.best_P    % Highest joint density
%   OPT.X_bar     % Sample mean
%   OPT.s         % Adaptive Metropolis scaling factor
%
%
%
%   Written by Khoa T. Tran, School of EE & CS
%                            University of Newcastle
%        		             Australia.
%
%   Reference: Andrieu, C., & Thoms, J. (2008). A tutorial on adaptive
%   MCMC. Statistics and Computing, 18(4), 343?373.
%
%   Definitions of additional terms in the codes
%   ---------------------------------
%   target    % Compute the log of the density function. To facilitate
%               Matlab coder translation, this function is to be hardcoded
%               by the user in a separate target.m file in the Matlab path.
%   delta     % Historical adaptation steps
%   acc_rate  % Metropolis acceptance rate

function OPT = am(target,OPT)
%% Initiation
coder.extrinsic('warning','cputime','num2str','save');
OPT.Dims        = length(OPT.X);        % Number of dimensions
OPT.best_X      = OPT.X;                % Record highest density state
OPT.P0          = target(OPT.X);        % Last density
OPT.best_P      = OPT.P0;               % Best seen density
OPT.eps         = -0.8;                 % Could be -0.5 to -1
OPT.burnin      = OPT.M;                % Initial burnin set to 100% sample
pcom            = 0;                    % Percentage complete count initialised to zero;
chk_stop        = 2000;                 % Periodically check to stop burn-in
delta           = zeros(OPT.burnin,2);  % Step change in adaptation
%% Take the Cholesky factor of V
[OPT.chol,errV]=chol(OPT.Sig);
if errV
    warning('Bad covariance matrix at start, perform reconditioning..')
    [V3,D] = eig(OPT.Sig);d=diag(D);d(d<=0)=10*realmin;
    OPT.Sig= V3*diag(d)*V3';
    L = tril(OPT.Sig,-1);OPT.Sig=diag(diag(OPT.Sig))+L+L';
    OPT.chol=chol(OPT.Sig);
end
%% Define optimal acceptance rate
switch OPT.Dims 
    case 1
        OPT.optimal=0.441;
    case 2
        OPT.optimal=0.352;
    case 3
        OPT.optimal=0.316;
    case 4
        OPT.optimal=0.285;
    case 5
        OPT.optimal=0.275;
    case 6
        OPT.optimal=0.273;
    case 7
        OPT.optimal=0.270;
    case 8
        OPT.optimal=0.268;
    case 9
        OPT.optimal=0.267;
    case 10
        OPT.optimal=0.266;
    case 11
        OPT.optimal=0.265;
    case 12
        OPT.optimal=0.264;
    otherwise
        OPT.optimal=0.255;
end
% Set initial scaling factor
OPT.s = 2.38^2/OPT.Dims;
%% Adaptive Metropolis Algorithm.
U = log(rand(OPT.burnin,1));        % Random acceptance numbers 
RD = randn(OPT.burnin,OPT.Dims);    % Random distubance proposal
prerun = 100*OPT.Dims;              % Generate 'prerun' samples before burnin
X_2 = zeros(prerun,OPT.Dims);       % Pre-burnin samples
OPT.avg_acc = 0;                    % Average acceptance rate
%% Pre-burnin run - generates 'prerun' samples with the initial covariance matrix
mark = cputime;mark_0 = cputime;% Used to keep track of elapsed time
for i = 1:prerun
    if OPT.dsp % Display completion status
        if mod(i,floor(OPT.burnin/20))==0
            pcom=pcom+5;
            fprintf('Complete = %d%% burn-in, Time since last update = %f s\n',pcom,cputime-mark); 
            if i>2
            remaining = (20 - i/floor(OPT.burnin/20))*(cputime-mark);
            hrs = floor(remaining/3600); remaining = rem(remaining,3600);
            mins = floor(remaining/60); 
            secs = floor(rem(remaining,60));
            fprintf('Predicted burn-in completion in %d hrs:%d mins:%d secs\n',hrs,mins,secs)
            fprintf('Acceptance rate %f %%\n',100*OPT.avg_acc);
            end;
        mark=cputime;
        end;
    end; 
    X_new = OPT.X + sqrt(OPT.s)*RD(i,:)*OPT.chol; % propose a new sample
    P_new = target(X_new);          % new density
    if P_new>OPT.best_P             % new density is highest 
        OPT.best_X=X_new;           % record best seen parameters
        OPT.best_P = P_new;
        OPT.X = X_new;              % accept X_new 
        OPT.P0 = P_new;
        acc_rate=1;                 % record acceptance rate
    else 
        rho = P_new - OPT.P0;
        acc_rate=exp(min(rho,0));   % record acceptance rate
        if U(i)<= rho
        OPT.X = X_new;              % accept X_new 
        OPT.P0 = P_new;
        end
    end
    X_2(i,:) = OPT.X;               % Store prerun samples
    %Update scale factor by Robbins-Monro approximation and record everage acceptance rate
    delta(i,1) = (i^OPT.eps)*(acc_rate-OPT.optimal);
    OPT.s = exp(log(OPT.s)+delta(i,1));
    OPT.avg_acc = OPT.avg_acc*(i-1)/i + acc_rate/i;
end
OPT.X_bar=mean(X_2);            % Sample mean
if OPT.cov == 1
    OPT.Sig=cov(X_2);             % Sample covariance matrix
    L = tril(OPT.Sig,-1);OPT.Sig=diag(diag(OPT.Sig))+L+L'; % Check for symmetry
    [OPT.chol,errV]=chol(OPT.Sig);
    if errV
        warning('Bad covariance matrix at prerun, perform reconditioning..')
        [V3,D] = eig(OPT.Sig);d=diag(D);d(d<=0)=10*realmin;
        OPT.Sig= V3*diag(d)*V3';
        L = tril(OPT.Sig,-1);OPT.Sig=diag(diag(OPT.Sig))+L+L'; % Check for symmetry again
        OPT.chol=chol(OPT.Sig);
    end
end
i=prerun+1;
%% Burnin run
while i <= OPT.burnin
    if OPT.dsp % Display completion status
        if mod(i,floor(OPT.burnin/20))==0
            pcom=pcom+5;
            fprintf('Complete = %d%% burn-in, Time since last update = %f s\n',pcom,cputime-mark); 
            if i>2
            remaining = (20 - i/floor(OPT.burnin/20))*(cputime-mark);
            hrs = floor(remaining/3600); remaining = rem(remaining,3600);
            mins = floor(remaining/60); 
            secs = floor(rem(remaining,60));
            fprintf('Predicted burn-in completion in %d hrs:%d mins:%d secs\n',hrs,mins,secs)
            fprintf('Acceptance rate %f %%\n',100*OPT.avg_acc);
            end;
        mark=cputime;
        end;
    end; 
    X_new = OPT.X + sqrt(OPT.s)*RD(i,:)*OPT.chol;     % propose a new sample
    P_new = target(X_new);              % new density
    if P_new>OPT.best_P                 % new density is highest 
        OPT.best_X=X_new;               % record best seen parameters
        OPT.best_P = P_new;
        OPT.X = X_new;                  % accept X_new
        OPT.P0 = P_new;
        acc_rate=1;                     % record acceptance rate
    else 
        rho = P_new - OPT.P0;
        acc_rate=exp(min(rho,0));       % record acceptance rate
        if U(i)<= rho
        OPT.X = X_new;                  % accept X_new 
        OPT.P0 = P_new;
        end
    end
    % Check stopping rule at every chk_stop counts
    if ~mod(i-prerun,chk_stop)
        if (log(mean(abs(delta(i-chk_stop+1:i-chk_stop/2,1)))) <= log(mean(abs(delta(i-chk_stop/2:i-1,1))))) &&...
           (log(mean(abs(delta(i-chk_stop+1:i-chk_stop/2,2)))) <= log(mean(abs(delta(i-chk_stop/2:i-1,2)))))
            OPT.burnin = i;
            delta = delta(1:OPT.burnin,:);   % Adaptation cost
        end
    end
    % Update scale
    gamma = i^OPT.eps;
    delta(i,1) = gamma*(acc_rate-OPT.optimal);
    OPT.s = exp(log(OPT.s)+delta(i,1));
    OPT.avg_acc = OPT.avg_acc*(i-1)/i + acc_rate/i;
    % Update Sigma
    if OPT.cov == 1
        delta_X=(OPT.X-OPT.X_bar); OPT.X_bar = OPT.X_bar+gamma*delta_X;
        Sdelta=gamma*(delta_X'*delta_X)-((i-1)^OPT.eps)*OPT.Sig;
        delta(i,2) = norm(Sdelta,1);
        OPT.Sig = OPT.Sig+Sdelta; L = tril(OPT.Sig,-1);OPT.Sig=diag(diag(OPT.Sig))+L+L'; % This help keeping V symmetrical
        [OPT.chol,errV]=chol(OPT.Sig);
        if errV
            warning('Bad covariance matrix at burnin, perform reconditioning..')
            [V3,D] = eig(OPT.Sig);d=diag(D);d(d<=0)=10*realmin;
            OPT.Sig= V3*diag(d)*V3';
            L = tril(OPT.Sig,-1);OPT.Sig=diag(diag(OPT.Sig))+L+L';
            OPT.chol=chol(OPT.Sig);
            i=i-1;
        end
    end
    i=i+1;
end
OPT.cputime.burnin = cputime-mark_0; % Record time taken to adapt
if OPT.dsp % If requested, give feedback on completion status
    if OPT.cov == 1
        figure;plot(1:OPT.burnin,delta(:,2));
        title('Adaptation process of the covariance matrix $$\mathbf{\Sigma}$$');set(gca,'YScale','log');
        xlabel('Number of iterations');ylabel('$$\Delta_m = \Vert\Sigma_m - \Sigma_{m-1}\Vert_1$$');
        axis tight;drawnow
    end
    
    figure;plot(1:OPT.burnin,abs(delta(:,1)));
    title('Adaptation process of the scaling factor $$s$$');set(gca,'YScale','log');
    xlabel('Number of iterations');ylabel('$$\Delta_m = \Vert s_m - s_{m-1}\Vert$$');
    axis tight;drawnow
    
    fprintf('Burn-in completed after %d samples\n',OPT.burnin);
    fprintf('Acceptance rate %f %%\n',100*OPT.avg_acc);
end
%% Final run
OPT.div = ceil(OPT.M/1e5);                  % Generate 'div' consecutive set of samples
OPT.mrun = ceil(OPT.M/OPT.div);             % with 'mrun' samples in each set
his = zeros(1.2*OPT.hgrd+1,OPT.Dims);   % Histogram grid points
OPT.X_bar = zeros(size(OPT.X));         % Sample mean vector
pcom = 0;                               % Completion percentage 
mark = cputime;mark_0 = mark;           % Used to keep track of elapsed time
OPT.TH = zeros(OPT.mrun,OPT.Dims);      % Final MCMC output
for m=1:OPT.div                             % Devide run into 'div' number of pieces to save RAM
    OPT.m = m;
    OPT.avg_acc = 0;                    % Average acceptance rate
%     keyboard
    OPT = am_kernel(target,OPT);
    %% Define grid points
    h_max = max(OPT.TH);    h_min = min(OPT.TH);    
    if m==1 % First grid
        range = 1.2;
        h_range = h_max-h_min;
        hist_grid = my_linspace(h_min-0.1*h_range,h_max+0.1*h_range,range*OPT.hgrd+1);
        mark_low = min(hist_grid,[],3); mark_high = max(hist_grid,[],3); 
    else    % Subsequence grid
        while any(h_max>mark_high) || any(h_min<mark_low)
            % Expand grid points
            range = range + 0.2;
            hist_grid = my_linspace(mark_low-0.1*h_range,mark_high+0.1*h_range,range*OPT.hgrd+1);
            mark_low = min(hist_grid,[],3); mark_high = max(hist_grid,[],3); 
            % Added zeros into old histogram
            his = [zeros(0.1*OPT.hgrd,OPT.Dims); his; zeros(0.1*OPT.hgrd,OPT.Dims)]; %#ok<AGROW>
        end
    end
    %% Compute raw sample histogram
    for k=1:OPT.Dims
        [memo,~] = hist(OPT.TH(:,k),hist_grid(1,k,:)); 
        his(:,k)= his(:,k)*(m-1)/m+memo'/m;
    end
    %% Compute sample mean
    OPT.X_bar = OPT.X_bar*(m-1)/m+mean(OPT.TH)/m;
    %% Save MCMC samples to hard disk and free RAM
    save([OPT.filename '_' num2str(m)],'OPT','-v7.3')
    if OPT.dsp  % Display completion status
        pcom=pcom+100/OPT.div;
        fprintf('Complete = %d%%, Time since last update = %f s\n',pcom,cputime-mark); 
        remaining = (OPT.div - m)*(cputime-mark);
        hrs = floor(remaining/3600); remaining = rem(remaining,3600);
        mins = floor(remaining/60); 
        secs = floor(rem(remaining,60));
        fprintf('Predicted completion in %d hrs:%d mins:%d secs\n',hrs,mins,secs)
        fprintf('Acceptance rate %f %%\n',100*OPT.avg_acc);
        mark=cputime;
    end
end
OPT.cputime.run = cputime-mark_0;         % Record time taken to run MCMC
% Normalise sample histograms in each dimension
for k=1:OPT.Dims
    h.p = his(:,k);h.x = hist_grid(1,k,:);h.x = h.x(:);
    h.p = h.p/(sum(h.p)*(h.x(2)-h.x(1)));
    OPT.hist.p(:,k) = h.p; OPT.hist.x(:,k) = h.x;
end
if OPT.dsp  % Display final results
    fprintf('Successfully generated %d samples in %d seconds\n',OPT.div*OPT.mrun, OPT.cputime.run);
end