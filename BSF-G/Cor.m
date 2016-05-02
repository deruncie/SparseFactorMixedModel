function [ r ] = Cor( X,Y )
%correlation matrix of correlations between rows of X and Y

x=size(X,1);
y=size(Y,1);

X=[X;Y];

covX=cov(X');
d=diag(sqrt(1./diag(covX)));
r=d *covX*d;

r=r(1:x,x+1:x+y);


end

