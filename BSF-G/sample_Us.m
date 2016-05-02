function [genetic_effects] = sample_Us(Factors,genetic_effects,svd_ZZ_Ainv)
%samples genetic effects (U) conditional on the factor scores F:
% F_i = U_i + E_i, E_i~N(0,s2*(h2*ZAZ + (1-h2)*I)) for each latent trait i
% U_i = zeros(r,1) if h2_i = 0
% it is assumed that s2 = 1 because this scaling factor is absorbed in
% Lambda
% svd_ZZ_Ainv has parameters to diagonalize a*Z_1*Z_1' + b*I for fast
% inversion:



Q = svd_ZZ_Ainv.Q;
s1 = svd_ZZ_Ainv.s1;
s2 = svd_ZZ_Ainv.s2;


global Z_1
k = Factors.k;
n = genetic_effects.n;
tau_e = 1./(1-Factors.h2);
tau_u = 1./Factors.h2;
b = Q'*Z_1*bsxfun(@times,Factors.scores,tau_e)';
z = randn(n,k);
for j=1:k
    if tau_e(j)==1
        genetic_effects.U(j,:) = zeros(1,n);
    elseif tau_e(j) == Inf
            genetic_effects.U(j,:) = Factors.scores(j,:);
    else
        d = s2*tau_u(j) + s1*tau_e(j);
        mlam = bsxfun(@times,b(:,j),1./d);
        genetic_effects.U(j,:) = (Q*(mlam + bsxfun(@times,z(:,j),1./sqrt(d))))';
    end
end
end


