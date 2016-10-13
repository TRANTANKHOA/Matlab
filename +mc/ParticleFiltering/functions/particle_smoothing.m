function Path = particle_smoothing( Path, theta, rho )
%PARTICLE_SMOOTHING Recalculating the weights of the particle filter result
%to reflect the backward smoothing information
T = length(Path);   
m = 1e4;
X1 = important_resampling(Path(end).state,Path(end).w,m);
w1 = 1/m*ones(m,1);
Path(end).state = X1;
Path(end).w = w1;
Path(end).ESS = 1;
formatSpec = 'Completed step %ith of %i steps\n';

for l=1:T-1
    X0 = important_resampling(Path(end-l).state,Path(end-l).w,m);
    M  = smoothing_matrix( X1,X0,theta );
    w1 = sum( bsxfun(@times,M,m*w1'./sum(M,1)), 2 );
    w1 = w1/sum(w1);
    Path(end-l).ESS = 1/sum(w1.^2)/m;
    %% Resampling if necessary
    if Path(end-l).ESS<=rho
        X0 = important_resampling( X0,w1,m );
        w1(:) = 1/m;
    end
    %     Path(end-l).ESS = 1;
    Path(end-l).w = w1;
    Path(end-l).state = X0;
    X1 = X0;
    clc
    fprintf(formatSpec,l,T)
end