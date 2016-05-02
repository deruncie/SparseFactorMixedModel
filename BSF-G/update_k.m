function  [F,Lambda,F_a,F_h2,Lambda_prec,Plam,delta,tauh] = update_k( F,Lambda,F_a,F_h2,Lambda_prec,Plam,delta,tauh,Z_W,Lambda_df,delta_2_shape,delta_2_rate,b0,b1,i,epsilon,prop )
%adapt the number of factors by dropping factors with only small loadings
%if they exist, or adding new factors sampled from the prior if all factors
%appear important. The rate of adaptation decreases through the chain,
%controlled by b0 and b1. Should work correctly over continuations of
%previously stopped chains.

[n,k] = size(F);
r = size(F_a,1);
p = size(Lambda,1);
gene_rows = 1:p;

prob = 1/exp(b0 + b1*i);                % probability of adapting
uu = rand;
lind = mean(abs(Lambda(gene_rows,:)) < epsilon);    % proportion of elements in each column less than eps in magnitude
vec = lind >=prop;num = sum(vec);       % number of redundant columns


if uu < prob && i>200
    if  i > 20 && num == 0 && all(lind < 0.995) && k < 2*p %add a column
        k=k+1;
        Lambda_prec(:,k) = gamrnd(Lambda_df/2,2/Lambda_df,[p,1]);
        delta(k,1) = gamrnd(delta_2_shape,1/delta_2_rate);
        tauh = cumprod(delta);
        Plam = bsxfun(@times,Lambda_prec,tauh');
        Lambda(:,k) = randn(p,1).*sqrt(1./Plam(:,k));
        F_h2(k) = rand;
        F_a(:,k) = randn(r,1)*sqrt(F_h2(k));
        F(:,k) = Z_W*F_a(:,k) + randn(n,1).*sqrt(1-F_h2(k));
    elseif num > 0      % drop redundant columns
        nonred = setdiff(1:k,find(vec)); % non-redundant loadings columns
        k = max(k - num,1);
        Lambda = Lambda(:,nonred);
        Lambda_prec = Lambda_prec(:,nonred);
        F = F(:,nonred);
        for red = setdiff(1:k-1,nonred)
            %combine deltas so that the shrinkage of kept columns doesnt
            %decrease after dropping redundant columns
            delta(red+1) = delta(red+1)*delta(red);
        end
        delta = delta(nonred);
        tauh = cumprod(delta);
        Plam = bsxfun(@times,Lambda_prec,tauh');
        F_h2 = F_h2(nonred);
        F_a = F_a(:,nonred);
    end
end


end

