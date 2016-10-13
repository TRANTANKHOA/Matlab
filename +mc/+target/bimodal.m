function log_p = bimodal(theta)
%LOG_DENSITY Give the log density of the point theta according to a bimodal
% target density. Each point is stored on a row.

% Some constants
Dims = size(theta,2);          % Length of theta
% Lehmer's positive definite matrix of order Dims
M = lehmer(Dims);
N = inv(M);
d = 10;
% Define 2 modes with multivariate Gaussian density
x1 = log_mvnpdf(theta,0*ones(1,Dims),N);   % 1st mode has zero mean and 
                                           % identity covariance
x2 = log_mvnpdf(theta,d*ones(1,Dims),M);   % 2nd mode 

% Calculate the combined density on the log scale plus a normalising
% constant (Z = Dims) in the end.
log_p = lwse_row([x1,x2],[0,1])+Dims;



