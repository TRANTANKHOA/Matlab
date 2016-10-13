function M = lehmer(Dims)
% Giving the Lehmer's positive definite matrix of order Dims
M = ones(Dims,1)*(1:Dims);
M = M./M';
M = tril(M) + tril(M,-1)';