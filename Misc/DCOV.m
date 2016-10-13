function [ CX ] = DCOV( X,px )
%Calculate covariance matrix CX for discrete 
%probability mass function pX of random vector X
%pX is a ROW vector and each ROW of X is a sample
EX = px*X;
nmeas = size(X,1);
XB = X - repmat(EX,nmeas,1); %centre to 0
PD = zeros(nmeas, nmeas);
PD(logical(eye(size(PD)))) = px; 
CX = XB'*PD*XB;
end

