U1 = dir('u1oversampled*.ext');
N1 = length(U1);
L2_1 = zeros(N1,1);
for k = 1:N1
    u1 = load(U1(k).name);
    n1 = length(u1)-1;
    h1 = 1.0/n1;
    x1 = 0 + h1*(0:n1)';
    uexact1 = -sin(4*pi*x1);
    L2_1(k) = sqrt(h1*sum((u1(1:n1/2+1) - uexact1(1:n1/2+1)).^2));
end

U2 = dir('u2oversampled*.ext');
N2 = length(U2);
L2_2 = zeros(N2,1);
for k = 1:N2
    u2 = load(U2(k).name);
    n2 = length(u2)-1;
    h2 = 1.0/n2;
    x2 = 0 + h2*(0:n2)';
    uexact2 = sin(2*pi*x2);
    L2_2(k) = sqrt(h2*sum((u2(n2/2+1:end) - uexact2(n2/2+1:end)).^2));
end

L2 = L2_1 + L2_2;
L2conv = log2(L2(1:end-1)./L2(2:end))


