function mass= simpson2array(array,grid)
%% Take a 1D row array and grid point to compute integral with Simpson's rule
array = array(:)';
grid  = grid(:)';
n = length(array);
if mod(n,2)==0
    array = [array 0];
    range = (max(grid(:))-min(grid(:)))*n/(n-1);
    n=n+1;
else
    range = (max(grid(:))-min(grid(:)));
end
weight = range/(3*(n-1));
k = 1:1:n-2; 
w=[1,2.^(mod(k,2)+ones(size(k))),1];  % Use Simpson's Rule   
mass = weight*sum(w.*array);
end