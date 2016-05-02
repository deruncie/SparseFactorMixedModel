function [Lambda] = sample_Lambda( Wtil,F,resid_W_prec, E_a_prec,Plam,invert_aI_bZAZ )
%Sample factor loadings Lambda while marginalizing over residual
%genetic effects: Y - Z_2W = F*Lambda' + E, vec(E)~N(0,kron(Psi_E,In) + kron(Psi_U, ZAZ^T))
%note: conditioning on F, but marginalizing over E_a.
%sampling is done separately by trait because each column of Lambda is
%independent in the conditional posterior
%note: invert_aI_bZAZ has parameters that diagonalize aI + bZAZ for fast
%inversion: inv(aI + bZAZ) = 1/b*U*diag(1./(s+a/b))*U'

p=size(resid_W_prec,1);
k=size(F,2);

U = invert_aI_bZAZ.U;
s = invert_aI_bZAZ.s;
FtU = F'*U;
UtY = U'*Wtil;

Zlams = randn(k,p);
Lambda = zeros(p,k);
for j = 1:p,
    FUDi = E_a_prec(j)*bsxfun(@times,FtU,1./(s' + E_a_prec(j)/resid_W_prec(j)));
    means = FUDi*UtY(:,j);
    Qlam = FUDi*FtU' + diag(Plam(j,:)); 
    Llam = chol(Qlam,'lower');
    vlam = Llam\means; mlam = Llam'\vlam; ylam = Llam'\Zlams(:,j);
    Lambda(j,:) = (ylam + mlam);
end


% save('sample_Lambda_data')
end
