function [mu,sigma,s,h,r] = ell_alp(alpha,theta,y)
    tau_2 = exp(theta);
    sigma = sqrt(tau_2./(1+tau_2));
    mu    = y.*sigma^2;
    s     = exp(-y.^2./(2.*...
        (tau_2+1)))./sqrt(tau_2+1);
    r = -0.5.*(alpha-mu).^2./sigma.^2;
    h = r - log(sigma) + log(s);
end
