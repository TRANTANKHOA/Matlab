function [X1,idx] = important_resampling( X0,w,N )
%IMPORTANT_RESAMPLING Perform a systematic resampling from the particles X0
%with weights w.
%
%Each particle is stored on a row. 
u = rand/N+(0:1/N:1-1/N)';
F = cumsum(w);
idx = zeros(N,1); 
tab = 1;
for i=1:N
    while F(tab)<=u(i)
        tab = tab+1;
    end
    idx(i) = tab;
end
X1 = X0(idx,:);
