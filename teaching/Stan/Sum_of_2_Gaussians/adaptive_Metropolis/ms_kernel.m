function OPT = ms_kernel(target,OPT)
% Random [0,1] numbers on log OPT.s
U = exprnd(1,OPT.mrun,1);       
% Uniform random vectors from ellipsoid
RD = randn(OPT.mrun,OPT.Dims);
RD = sqrt(OPT.s)*bsxfun(@times,...
                    bsxfun(@rdivide,...
                        RD,sqrt(sum(RD.^2,2))...
                    ),rand(OPT.mrun,1).^(1/OPT.Dims)...
                 )*OPT.chol;    
% Initial V
Sigma = OPT.s*OPT.Sig;
OPT.V     = sqrt(OPT.s)*randn(1,OPT.Dims)*OPT.chol;   
OPT.PV    = log_mvnpdf(OPT.V,zeros(1,OPT.Dims),Sigma);
OPT.Jcb   = log(det(Sigma));
for i=1:OPT.mrun
    %% Update h
    log_h = OPT.P0 + OPT.PV - U(i);
    %% Update V
    rho    = -2*(log_h - OPT.P0) - OPT.Dims*log(2*pi) - OPT.Jcb;
    OPT.V  = RD(i,:).*sqrt(rho);
    OPT.PV = log_mvnpdf(OPT.V,zeros(1,OPT.Dims),Sigma);
    %% Update X
    X_new = OPT.X + OPT.V;    % propose a new sample
    P_new = target(X_new);          % new density
    if (P_new + OPT.PV)> log_h                 % new density is highest 
        OPT.X = X_new; OPT.P0 = P_new;
        acc_rate = 1;
    else
        acc_rate = 0;
    end
    OPT.TH(i,:) = OPT.X;
    OPT.avg_acc = OPT.avg_acc*(i-1)/i + acc_rate/i;
end