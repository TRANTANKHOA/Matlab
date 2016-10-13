function p=lnNpr(a,b)
% computes ln(P(a<Z<b))
% where Z~N(0,1) very accurately for any 'a', 'b'
p=zeros(size(a));
% case b>a>0
I=a>0;
if any(I)
    pa=lnPhi(a(I)); % log of upper tail
    pb=lnPhi(b(I));
    p(I)=pa+log1p(-exp(pb-pa));
end
% case a<b<0
idx=b<0;
if any(idx)
    pa=lnPhi(-a(idx)); % log of lower tail
    pb=lnPhi(-b(idx));
    p(idx)=pb+log1p(-exp(pa-pb));
end
% case a<0<b
I=(~I)&(~idx);
if any(I)
    pa=erfc(-a(I)/sqrt(2))/2; % lower tail
    pb=erfc(b(I)/sqrt(2))/2;  % upper tail
    p(I)=log1p(-pa-pb);
end