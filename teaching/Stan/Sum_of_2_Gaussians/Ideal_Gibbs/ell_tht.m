function h = ell_tht(alpha,theta,y)
    tau_2 = exp(theta);
    sigma = sqrt(tau_2./(1+tau_2));
    mu    = y.*sigma^2;
    s     = exp(-y.^2./(2.*...
        (tau_2+1)))./sqrt(tau_2+1);
    h = -0.5.*(alpha-mu).^2./sigma.^2 - log(sigma) + log(s);
end

