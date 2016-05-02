function  [ Factors,genetic_effects] = update_k( Factors,genetic_effects,b0,b1,i,epsilon,prop )
%adapt the number of factors by dropping factors with only small loadings
%if they exist, or adding new factors sampled from the prior if all factors
%appear important. The rate of adaptation decreases through the chain,
%controlled by b0 and b1

global Z_1

df = Factors.df;
ad2 = Factors.ad2;
bd2 = Factors.bd2;
p = Factors.p;
k = Factors.k;
r = Factors.r;
gene_rows = 1:p;
Lambda = Factors.Lambda;

prob = 1/exp(b0 + b1*i);                % probability of adapting
uu = rand;
lind = mean(abs(Lambda(gene_rows,:)) < epsilon);    % proportion of elements in each column less than eps in magnitude
vec = lind >=prop;num = sum(vec);       % number of redundant columns

Factors.num=num;
Factors.no_f(i)=k-num;

if uu < prob && i>200
    if  i > 20 && num == 0 && all(lind < 0.995) && k < 2*p %add a column
        k=k+1;
        Factors.k = k;
        Factors.psijh(:,k) = gamrnd(df/2,2/df,[p,1]);
        Factors.delta(k,1) = gamrnd(ad2,1/bd2);
        Factors.tauh = cumprod(Factors.delta);
        Factors.Plam = bsxfun(@times,Factors.psijh,Factors.tauh');
        Factors.Lambda(:,k) = randn(p,1).*sqrt(1./Factors.Plam(:,k));
        Factors.h2(k) = rand;
        Factors.accept_h2_proposal(k) = 0;
        genetic_effects.U(k,:) = randn(1,genetic_effects.n);
        Factors.scores(k,:) = genetic_effects.U(k,:)*Z_1 + randn(1,r).*sqrt(1-Factors.h2(k));
    elseif num > 0      % drop redundant columns
        nonred = setdiff(1:k,find(vec)); % non-redundant loadings columns
        k = max(k - num,1);
        Factors.k = k;
        Factors.Lambda = Lambda(:,nonred);
        Factors.psijh = Factors.psijh(:,nonred);
        Factors.scores = Factors.scores(nonred,:);
        for red = setdiff(1:k-1,nonred)
            %combine deltas so that the shrinkage of kept columns doesnt
            %decrease after dropping redundant columns
            Factors.delta(red+1) = Factors.delta(red+1)*Factors.delta(red);
        end
        Factors.delta = Factors.delta(nonred);
        Factors.tauh = cumprod(Factors.delta);
        Factors.Plam = bsxfun(@times,Factors.psijh,Factors.tauh');
        Factors.h2 = Factors.h2(nonred);
        genetic_effects.U = genetic_effects.U(nonred,:);
    end
end
Factors.nofout(i+1)=k;


end

