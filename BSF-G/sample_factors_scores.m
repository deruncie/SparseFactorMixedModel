function [ F ] = sample_factors_scores( W_tilde, Z_W,Lambda,resid_W_prec,F_a,F_h2 )
%Sample factor scores given factor loadings (F_a), factor heritabilities (F_h2) and
%phenotype residuals
    Lmsg = bsxfun(@times,Lambda,resid_W_prec);
    tau_e = 1./(1-F_h2);
    S=chol(Lambda'*Lmsg+diag(tau_e),'lower');
    Meta = ((W_tilde*Lmsg + bsxfun(@times,Z_W*F_a,tau_e'))/S')/S;
    F = Meta + randn(size(Meta))/S;   
end

