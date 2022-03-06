
load('M')
load('P')
load('R')
load('B')
load('C')
Bnew = zeros(36);
for i = 1:36
    for j = 1:36
        Bnew(i,j) = M(i,j)*exp(R(i)+C(j));
    end
end
 
I = eye(36);
PERM = I(:,P);

Bnew = Bnew*PERM;
