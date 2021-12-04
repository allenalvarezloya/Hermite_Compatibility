load('Mu_inner.ext');
s = svd(Mu_inner);
[U,S,V] = svd(Mu_inner);
surf(U(:,end)*V(:,end)')
v = ones(length(Mu_inner),1);
RHS = Mu_inner*v;
conditionNumber = cond(Mu_inner)
svdSol = 0*Mu_inner;
for k = 1:1
    svdSol = V(:,1:k)*S(1:k,1:k).^(-1)*U(:,1:k)'*v;
end
    