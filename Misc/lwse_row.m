function s = lwse_row(x,w)
% Compute log(sum(exp(x).*w)) while avoiding numerical underflow.
% Each points of x is on a column and the weights w is on a row.
x = bsxfun(@plus,x,log(w)); % weighting the elements of each row of x with
                            % corresponding row in w
% Calculating the log(sum(exp(x).*w)) for rows with infinite maximum value
s = max(x,[],2);         % find the maximum value in each row
flags = isfinite(s);     % marking the finite elements in xmax

% Calculating the log(sum(exp(x).*w)) for rows with finite maximum value
x(flags,:) = bsxfun(@minus,x(flags,:),s(flags));
s(flags) = s(flags) + log(sum(exp(x(flags,:)),2));
