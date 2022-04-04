M = dir('M*.ext');
N = length(M);
MAX = 0;
for k = 1:N
    m = load(M(k).name);
    rank(m);
    sM = max(MAX,cond(m));
end


B = dir('B*.ext');
N = length(B);
MAX = 0;
for k = 1:N
    b = load(B(k).name);
    max(max(abs(b)));
    sB = max(MAX,cond(b));
end
sM 
sB