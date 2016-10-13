function [ M ] = smoothing_matrix( X1,X0,theta )
%SMOOTHING_MATRIX Calculate the matrix of conditional probablities
%M(i,j) = p(x1_j|x0_i)
%
%The specifics of this function depend on the particular dynamical system,
%i.e. this function is user-defined.

SIGMA = 10^theta;   
M  = -0.5 * ( bsxfun( @rdivide, bsxfun(@minus,X1',0.7*X0).^2 , SIGMA) + log(SIGMA) + log(2*pi));
        
end

