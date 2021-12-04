SingularVals = dir('Singular*.ext');
s = load(SingularVals(1).name);
S = zeros(length(s),length(SingularVals));
for k = 2:length(SingularVals)
    S(:,k) = load(SingularVals(k).name);
end
surf(S)