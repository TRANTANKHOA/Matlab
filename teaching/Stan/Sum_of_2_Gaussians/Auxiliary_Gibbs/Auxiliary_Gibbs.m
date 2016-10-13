% This script implement the T-factors Auxiliary Gibbs/Slice sampling 
% approach for the model
% 
% $$ \alpha_t \sim N(0,1); t=1,2,...T $$
% 
% $$ y_t \sim N(\alpha_t,\exp(\theta)) $$
% 
%% Initialisation 
% M       = 1e5;  % # of Markov samples
P       = log_prior(theta);
theta_M = zeros(M,1);   % each row is one observation/sample
alpha_M = zeros(M,T);
w       = 1;    % initial width for slice sampling
%% Gibbs sampling + Slice sampling
for k=1:M
    %% Computing product of two Gaussians
    [mu_T,sigma,s_T,h_T,~] = ell_alp(alpha_T,theta,y);
    %% Sampling h_t, t=1,2,...T
    h       = h_T - exprnd(1,1,T);
    %% Sampling alpha_t, t=1,2,...T
    alpha_T = mu_T + (rand(1,T).*2-1).*...
        sigma.*sqrt(-2.*log(sigma.*exp(h)./s_T));
    %% Slice sampling on psi(theta) with condition ell(theta)>h
    ell = @(x) ell_tht(alpha_T,x,y);
%     if any( ell(theta) < h )
%         keyboard % for debugging only
%     end
    h_tht = P-exprnd(1,1,1);   % Draw a random height for each sample on the log scale.
    %% Generate a random initial interval [a,b] for each sample
    b = rand;    a = b - 1;
    %% Expanding the upper bound of the inteval [a,b]
    flag = true;                  % Marking the status of expansion. 
                                        % The expansion will continue for one sample if a 
                                        % flag is true (1) and stop if otherwise (0).
    n_w = 100;
    n_b = randi(n_w+1,1,1)-1;           % Maximum number of expansions of the upper bound
    i=0;                                % Reset the counter for the number of expansions
    while flag
        x = theta+w*b; 
        flag = all([log_prior(x), ell(x), n_b] > [h_tht, h, i]);                 % is density >= h?
        b = b+flag;                               % b = b+1
        i = i+1;
    end
%     clc
%     formatSpec = 'Expansions on upper bound is %d. k = %d.\n';
%     fprintf(formatSpec,[i-1,k])
    %% Expanding the lower bound of the inteval [a,b] (similar as above)
    flag = true;
    n_a = n_w - n_b;                    % Maximum number of expansions of the lower bound
    i=0;    
    while flag
        x = theta+w*a;
        flag = all([log_prior(x), ell(x), n_a] > [h_tht, h, i]);                 % is density >= h?
        a = a-flag;
        i = i+1;
    end
%     formatSpec = 'Expansions on lower bound is %d. k = %d.\n';
%     fprintf(formatSpec,[i-1,k])
    %% Proposing and shrinking bounds
    flag = true;
    while flag
        innv = (b-a).*rand+a;
        x    = theta+w*innv; 
        P    = log_prior(x);
        flag = any([P, ell(x) ] < [h_tht, h ]);
        if innv>=0
            b = innv;
        else
            a = innv;
        end
    end
    w = b-a;
    theta = x;  
    theta_M(k) = theta;
    alpha_M(k,:) = alpha_T;
end
time_Gibbs = toc
%% Plotting theta
window = 'Auxiliary Gibbs Sampling';
f = plot_time_series(theta_M,'$\upsilon$',500,window);
subplot(2,2,1)
set(gca,'Ylim',YL_theta)
subplot(2,2,2)
set(gca,'Ylim',YL_theta)
legend('off')

    %% Printing
    w = 0.5; h = 0.35; w=round(w*1050); h=round(h*800); r=300; sz = [-1000 1000 w h];
    drawnow; undock(f); set(f,'Units','pixels','Position',sz);
    filename = '/Users/Eric/Dropbox/PhD/Reports/Literature Review/Pics/sum_2_Gauss_theta_Gibbs.pdf';
    export_fig(filename,f,'-a1',['-r',num2str(r)]);%, '-grey','-nocrop'
    saveas(f,[filename(1:end-4) '.fig'])
%% Plotting alpha_t
f = plot_time_series(alpha_M(:,t),['$x_{' num2str(t) '}$' ],500,window);
subplot(2,2,1)
set(gca,'Ylim',YL_alpha)
subplot(2,2,2)
set(gca,'Ylim',YL_alpha)
legend('off')

    %% Printing
    w = 0.5; h = 0.35; w=round(w*1050); h=round(h*800); r=300; sz = [-1000 1000 w h];
    drawnow; undock(f); set(f,'Units','pixels','Position',sz);
    filename = '/Users/Eric/Dropbox/PhD/Reports/Literature Review/Pics/sum_2_Gauss_alpha_Gibbs.pdf';
    export_fig(filename,f,'-a1',['-r',num2str(r)]);%, '-grey','-nocrop'
    saveas(f,[filename(1:end-4) '.fig'])
%% Plotting whole trajectory of alpha_T
plot_alpha_T
set(gca,'Ylim',YL_alpha_M)
%% Plotting generated data y
Y = alpha_M + bsxfun(@times,randn(M,T),sqrt(exp(theta_M)));
plot_y_T
legend('Location','northeastoutside')
set(gca,'Ylim',YL_y)
    %% Printing
    axis tight
    w = 1; h = 0.25; w=round(w*1050); h=round(h*800); r=300; sz = [-1000 1000 w h];
    drawnow; undock(f); set(f,'Units','pixels','Position',sz);
    filename = '/Users/Eric/Dropbox/PhD/Reports/Literature Review/Pics/sum_2_Gauss_y_Gibbs.pdf';
    export_fig(filename,f,'-a1',['-r',num2str(r)]);%, '-grey','-nocrop'
    saveas(f,[filename(1:end-4) '.fig'])
