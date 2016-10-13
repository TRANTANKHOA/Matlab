function s = lse(x)
%LSE Compute log(sum(exp(x))) while avoiding numerical underflow.
% The vector x here is the importance weights on the logarithmic scale.
% This function is required for numerical stability. The vector x can only
% be a column, but not a matrix.

x = x(:);                           % Reformat the input into column vector
xmax = max(x);                      % Find maximum value in x
if isfinite(xmax)                   % Check if maximum value is finite
    x = x-xmax;                     % Normalise the x vector on the log scale
    s = xmax + log(sum(exp(x)));    % Calculate the log(sum(exp(x))) function
else
    s = xmax; 
end