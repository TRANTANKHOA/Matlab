function [V,C,err] = make_symm(V)
%Take an assymmetric matrix and make it symmetric. Check if it is positive
%definite and return the Cholesky factor if it is. Else, return a diagonal
%matrix with the maximum eigenvalues
L       = tril(V,-1);
V       = diag(diag(V))+L+L';
[C,err] = chol(V);
if err
    warning('Bad covariance matrix, performed reconditioning.')
    [V,D]   = eig(V);d=diag(D);d(d<=0)=10*realmin;
    V       = real(V*diag(d)*V');   
    C       = chol(V);
end


