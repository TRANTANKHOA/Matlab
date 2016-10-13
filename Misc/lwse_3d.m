function s = lwse_3d(x,w,dim)
% Compute log(sum(exp(x).*w)) while avoiding numerical underflow.
% Note that the rows of x are weighted by the rows of w

% x=cat(3, posterior(bis_flg,:,1), prior(bis_flg,:,1));w=cat(3, lambda_1(bis_flg), 1-lambda_1(bis_flg))
if nargin <3
    dim = 2;
end

x = bsxfun(@plus,x,log(w));
xmax = max(x,[],dim);
flags = isfinite(xmax);
shape = ones(1,dim);shape(dim) = size(w,dim);
idx = repmat(flags,shape);
xmax2= repmat(xmax,shape);
x(idx) = x(idx)-xmax2(idx);
s = log(sum(exp(x),dim));
s(flags) = s(flags) + xmax(flags);
s(~flags) = xmax(~flags);

