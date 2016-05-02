function [Factors] = sample_h2s_discrete(Factors,eig_ZAZ)
%sample factor heritibilties from a discrete set on [0,1)
%prior places 50% of the weight at h2=0
%samples conditional on F, marginalizes over U.

Ur = eig_ZAZ.vectors;
eta = eig_ZAZ.values;

r=Factors.r;
k = Factors.k;
s=Factors.h2_divisions;

log_ps = zeros(k,s);
std_scores_b = Factors.scores*Ur;
for i=1:s,
    h2 = (i-1)/(s);
    std_scores = Factors.scores;
    if h2 > 0,
        std_scores = 1/sqrt(h2)*bsxfun(@times,std_scores_b,1./sqrt(eta+(1-h2)/h2)');
        det = sum(log((eta+(1-h2)/h2)*h2)/2);
    else       
        det = 0;
    end
    log_ps(:,i) = sum(log(normpdf(std_scores)),2)- det;
    if i==1,
        log_ps = log_ps + log(s-1);
    end
end
for j=1:k
    norm_factor = max(log_ps(j,:))+log(sum(exp(log_ps(j,:)-max(log_ps(j,:)))));
    ps_j = exp(log_ps(j,:) - norm_factor);
    log_ps(j,:) = ps_j;
    Factors.h2(j) = sum(rand>cumsum(ps_j))/(s);
end

end