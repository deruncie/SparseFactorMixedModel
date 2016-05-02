
function [F_h2] =  sample_h2s_discrete(F,h2_divisions,h2_priors,invert_aI_bZAZ)
%sample factor heritibilties from a discrete set on [0,1)
%prior places 50% of the weight at h2=0
%samples conditional on F, marginalizes over F_a.
%uses invert_aI_bZAZ.U and invert_aI_bZAZ.s to not have to invert aI + bZAZ
%each iteration.

U = invert_aI_bZAZ.U;
s = invert_aI_bZAZ.s;

k = size(F,2);
F_h2 = zeros(k,1);

log_ps = zeros(k,h2_divisions);
std_scores_b = F'*U;
for i=1:h2_divisions,
    h2 = (i-1)/(h2_divisions);
    std_scores = F;
    if h2 > 0,
        std_scores = 1/sqrt(h2)*bsxfun(@times,std_scores_b,1./sqrt(s+(1-h2)/h2)')';
        det = sum(log((s+(1-h2)/h2)*h2)/2);
    else       
        det = 0;
    end
    log_ps(:,i) = sum(log(normpdf(std_scores)),1)'- det + log(h2_priors(i));
end
for j=1:k
    norm_factor = max(log_ps(j,:))+log(sum(exp(log_ps(j,:)-max(log_ps(j,:)))));
    ps_j = exp(log_ps(j,:) - norm_factor);
    log_ps(j,:) = ps_j;
    F_h2(j) = sum(rand>cumsum(ps_j))/(h2_divisions);
end


% save('sample_h2s_discrete_data')
end