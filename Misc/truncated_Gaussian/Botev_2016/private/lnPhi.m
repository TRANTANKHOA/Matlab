function p=lnPhi(x)
% computes logarithm of  tail of Z~N(0,1) mitigating
% numerical roundoff errors;
p=-0.5*x.^2-log(2)+reallog(erfcx(x/sqrt(2)));