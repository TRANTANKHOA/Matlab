function x=Phinv(p,l,u)
% computes with precision the quantile function
% of the standard normal distribution,
% truncated to the interval [l,u], using erfcinv.
I=u<0;l(I)=-l(I);u(I)=-u(I); % use symmetry of normal
pu=erfc(u/sqrt(2))/2;pl=erfc(l/sqrt(2))/2;
x=sqrt(2)*erfcinv(2*(pl+(pu-pl).*p));
x(I)=-x(I); % adjust sign due to symmetry