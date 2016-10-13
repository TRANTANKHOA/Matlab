function xy = row_cov(X)
% Compute variance of each row in X, return a column vector.
m = size(X,2);
xc = bsxfun(@minus,X,sum(X,2)/m);  % Remove mean
xy = sum(xc.^2,2) / (m-1);

