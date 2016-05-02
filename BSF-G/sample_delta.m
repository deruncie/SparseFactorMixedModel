function [ delta,tauh ] = sample_delta( Factors,Lambda2_resid )
%sample delta and tauh parameters that control the magnitudes of higher
%index factor loadings.

ad1 = Factors.ad1;
ad2 = Factors.ad2;
bd1 = Factors.bd1;
bd2 = Factors.bd2;
k = Factors.k;
delta = Factors.delta;
tauh = Factors.tauh;
psijh = Factors.psijh;

mat = bsxfun(@times,psijh,Lambda2_resid);
n_genes = size(mat,1);
ad = ad1 + 0.5*n_genes*k; bd = bd1 + 0.5*(1/delta(1))*sum(tauh.*sum(mat)');
delta(1) = gamrnd(ad,1/bd);
tauh = cumprod(delta);

for h = 2:k
    ad = ad2 + 0.5*n_genes*(k-h+1); bd = bd2 + 0.5*(1/delta(h))*sum(tauh(h:end).*sum(mat(:,h:end))');
    delta(h,1) = gamrnd(ad,1/bd);
    tauh = cumprod(delta);
end


end

