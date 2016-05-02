function corX = CovToCor(X)
%Normalize a covariance matrix into a correlation matrix
    corX=diag(1./sqrt(diag(X)))*X*diag(1./sqrt(diag(X)));
end

