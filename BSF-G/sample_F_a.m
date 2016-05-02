function [F_a] = sample_F_a(F,Z_W,F_h2,invert_aZZt_Ainv)
%samples genetic effects on factors (F_a) conditional on the factor scores F:
% F_i = F_{a_i} + E_i, E_i~N(0,s2*(1-h2)*I) for each latent trait i
% U_i = zeros(r,1) if h2_i = 0
% it is assumed that s2 = 1 because this scaling factor is absorbed in
% Lambda
% invert_aZZt_Ainv has parameters to diagonalize a*Z_1*Z_1' + b*I for fast
% inversion:



U = invert_aZZt_Ainv.U;
s1 = invert_aZZt_Ainv.s1;
s2 = invert_aZZt_Ainv.s2;

k = size(F,2);
n = size(Z_W,2);
tau_e = 1./(1-F_h2);
tau_u = 1./F_h2;
b = U'*Z_W'*bsxfun(@times,F,tau_e');
z = randn(n,k);
F_a = zeros(n,k);
for j=1:k
    if tau_e(j)==1
        F_a(:,j) = zeros(n,1);
    elseif tau_e(j) == Inf
        F_a(:,j) = F(:,j);
    else
        d = s2*tau_u(j) + s1*tau_e(j);
        mlam = bsxfun(@times,b(:,j),1./d);
        F_a(:,j) = U*(mlam + bsxfun(@times,z(:,j),1./sqrt(d)));
    end
end
end


