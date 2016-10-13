function s = lwse(x,w)
% Compute log(sum(exp(x).*w)) while avoiding numerical underflow.
% Each points of x is on a row and the weights w is on a column.
x = bsxfun(@plus,x,log(w)); % weighting the elements of each row of x with
                            % corresponding weight in w
% Filling rows with infinite maximum value with infinity
s = max(x);              % find the maximum value in each column
flags = isfinite(s);     % marking the finite elements in xmax

% Calculating the log(sum(exp(x).*w)) for rows with finite maximum value
x(:,flags) = bsxfun(@minus,x(:,flags),s(flags));
s(flags) = s(flags) + log(sum(exp(x(:,flags))));
