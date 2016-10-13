function xy = col_cov(X)
% Compute variance of each column in X, return a row vector.
m = size(X,1);
xc = bsxfun(@minus,X,sum(X,1)/m);  % Remove mean
xy = sum(xc.^2,1) / (m-1);

