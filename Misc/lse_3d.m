function s = lse_3d(data,dim)
% Compute log(sum(exp(x))) while avoiding numerical underflow.
% data = log_omega;
shape = ones(1,dim);shape(dim) = size(data,dim);
s = max(data,[],dim);
flag = isfinite(s);
flags = repmat(flag,shape);
ss = repmat(s,shape);
data(flags) = data(flags)-ss(flags);
lse_term = log(sum(exp(data),dim));
s(flag) = s(flag) + lse_term(flag);
