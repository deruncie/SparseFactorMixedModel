function [location_sample] = sample_means( Y_tilde, ...
            resid_Y_prec, E_a_prec, invert_aPXA_bDesignDesignT )
%when used to sample [B;E_a]:
% W - F*Lambda' = X*B + Z_1*E_a + E, vec(E)~N(0,kron(Psi_E,In)). 
% Note: conditioning on F, Lambda and W.
%The vector [b_j;E_{a_j}] is sampled simultaneously. Each trait is sampled separately because their
%conditional posteriors factor into independent MVNs.
%note:invert_aPXA_bDesignDesignT has parameters to diagonalize mixed model equations for fast inversion: 
%inv(a*blkdiag(fixed_effects_prec*eye(b),Ainv) + b*[X Z_1]'[X Z_1]) = U*diag(1./(a.*s1+b.*s2))*U'
%Design_U = [X Z_1]*U, which doesn't change each iteration. 
%

U = invert_aPXA_bDesignDesignT.U;
s1 = invert_aPXA_bDesignDesignT.s1;
s2 = invert_aPXA_bDesignDesignT.s2;
Design_U = invert_aPXA_bDesignDesignT.Design_U;

p = size(resid_Y_prec,1);
n = size(Design_U,2);

means = bsxfun(@times,Design_U'*Y_tilde,resid_Y_prec');
location_sample = zeros(n,p);
Zlams = randn(n,p);
for j = 1:p,
    d = s1*E_a_prec(j) + s2*resid_Y_prec(j);
    mlam = bsxfun(@times,means(:,j),1./d);
    location_sample(:,j) = U*(mlam + bsxfun(@times,Zlams(:,j),1./sqrt(d)));
end
location_sample=location_sample';

% save('sample_means_data')
end
