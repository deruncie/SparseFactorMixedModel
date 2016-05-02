function [ Factors ] = sample_factors_scores( Ytil, Factors,resid,genetic_effects )
%Sample factor scores given factor loadings, U, factor heritabilities and
%phenotype residuals
global Z_1
k = Factors.k;
n = Factors.n;
Lambda = Factors.Lambda;
Lmsg = bsxfun(@times,Lambda,resid.ps);
tau_e = 1./(1-Factors.h2);
S=chol(Lambda'*Lmsg+diag(tau_e),'lower');
Meta = S'\(S\(Lmsg'*Ytil + bsxfun(@times,genetic_effects.U*Z_1,tau_e)));
Factors.scores = Meta + S'\randn(k,n);   
end

