M = dir('M*.ext');
B = dir('B*.ext');
N = length(M);
MAX = 0;
for k = 1:N
    m = load(M(k).name);
    rank(m);
    sM = max(MAX,cond(m));
    b = load(B(k).name);
    max(max(abs(b)));
    sB = max(MAX,cond(b));
end
sM 
sB