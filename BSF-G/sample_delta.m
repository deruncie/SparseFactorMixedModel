function [ delta,tauh ] = sample_delta( delta,tauh,Lambda_prec,...
            delta_1_shape,delta_1_rate,delta_2_shape,delta_2_rate,Lambda2 )
%sample delta and tauh parameters that control the magnitudes of higher
%index factor loadings.

k = size(tauh,1);
mat = bsxfun(@times,Lambda_prec,Lambda2);
n_genes = size(mat,1);
shape = delta_1_shape + 0.5*n_genes*k; rate = delta_1_rate + 0.5*(1/delta(1))*sum(tauh.*sum(mat)');
delta(1) = gamrnd(shape,1/rate);
tauh = cumprod(delta);

for h = 2:k
    shape = delta_2_shape + 0.5*n_genes*(k-h+1); rate = delta_2_rate + 0.5*(1/delta(h))*sum(tauh(h:end).*sum(mat(:,h:end))');
    delta(h,1) = gamrnd(shape,1./rate);
    tauh = cumprod(delta);
end


% save('sample_delta_data')
end

