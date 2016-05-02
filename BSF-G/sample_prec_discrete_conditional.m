
function [Prec] =  sample_prec_discrete_conditional(Y,h2_divisions,h2_priors,invert_aI_bZAZ,res_prec)
%sample factor heritibilties conditional on a given residual precision
%(res_precision)
%prior given as relative weights on each of h2_divisions points. Doesn't
%have to be normalized
%samples conditional on F, marginalizes over F_a.
%uses invert_aI_bZAZ.U and invert_aI_bZAZ.s to not have to invert aI + bZAZ
%each iteration.
%ident_prec is a in the above equation.

U = invert_aI_bZAZ.U;
s = invert_aI_bZAZ.s;

n = size(Y,2);
Trait_h2 = zeros(n,1);

log_ps = zeros(n,h2_divisions);
std_scores_b = Y'*U;
for i=1:h2_divisions,
    h2 = (i-1)/(h2_divisions);
    std_scores = Y;
    if h2 > 0,
        std_scores = bsxfun(@times,bsxfun(@times,std_scores_b,1./sqrt(s+(1-h2)/h2)')',1./sqrt(h2./(res_prec'*(1-h2))));
        det = sum(log((s+(1-h2)/h2)*(h2./(res_prec'*(1-h2))))/2);
    else       
        det = 0;
    end
    log_ps(:,i) = sum(log(normpdf(std_scores)),1)'- det' + log(h2_priors(i));
end
for j=1:n
    norm_factor = max(log_ps(j,:))+log(sum(exp(log_ps(j,:)-max(log_ps(j,:)))));
    ps_j = exp(log_ps(j,:) - norm_factor);
    log_ps(j,:) = ps_j;
    Trait_h2(j) = sum(rand>cumsum(ps_j))/(h2_divisions);
end

Prec = (res_prec.*(1-Trait_h2))./Trait_h2;


% save('sample_prec_discrete_conditional_data')

end