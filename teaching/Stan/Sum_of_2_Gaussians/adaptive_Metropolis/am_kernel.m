function OPT = am_kernel(target,OPT)

U = -exprnd(1,OPT.mrun,1);              % Random acceptance numbers 
RD = randn(OPT.mrun,OPT.Dims)*OPT.chol; % Random distubance proposal
for i=1:OPT.mrun
    X_new = OPT.X + sqrt(OPT.s)*RD(i,:);% propose a new sample
    P_new = target(X_new);              % new density
    if P_new>OPT.best_P                 % new density is highest 
        OPT.best_X =X_new;              % record best seen parameters
        OPT.best_P = P_new;
        OPT.X      = X_new;             % accept X_new 
        OPT.P0 = P_new;
        acc_rate = 1;                     % record acceptance rate
    else 
        rho = P_new - OPT.P0;
        acc_rate=exp(min(rho,0));       % record acceptance rate
        if U(i)<= rho
        OPT.X = X_new;                  % accept X_new 
        OPT.P0 = P_new;
        end
    end
    OPT.TH(i,:) = OPT.X;
    %Update scale
    idx   = (i+(OPT.m-1)*OPT.mrun+OPT.burnin);    % Global index
    gamma = idx^OPT.eps;
    OPT.s = exp(log(OPT.s) + (acc_rate-OPT.optimal));
    OPT.avg_acc = OPT.avg_acc*(i-1)/i + acc_rate/i;
    %Update Sigma
    delta_X=OPT.X-OPT.X_bar;    OPT.X_bar = OPT.X_bar+gamma*delta_X;
    Sdelta=gamma*(delta_X'*delta_X)-((idx-1)^OPT.eps)*OPT.Sig;
    OPT.Sig = OPT.Sig+Sdelta;   L = tril(OPT.Sig,-1);OPT.Sig=diag(diag(OPT.Sig))+L+L';
    [OPT.chol,errV]=chol(OPT.Sig);
    if errV
        warning('Bad covariance matrix at generation state, perform reconditioning..')
        [V3,D] = eig(OPT.Sig);d=diag(D);d(d<=0)=10*realmin;
        OPT.Sig= V3*diag(d)*V3';
        L = tril(OPT.Sig,-1);OPT.Sig=diag(diag(OPT.Sig))+L+L';
        OPT.chol=chol(OPT.Sig);
    end
end