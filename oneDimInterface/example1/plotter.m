clear, clc
figure(1)
S1 = dir('solution1*.ext');
N = length(S1);
s1 = load(S1(1).name);
n1 = length(s1) - 1;
h1 = 1.0/n1;
x1 = 0.0 + (0:n1)'*h1;

for i = 1:N
    s1 = load(S1(i).name);
    plot(x1,s1,'LineWidth',2)
    axis([0 1 -1 1])
    TITLE = sprintf('%i',i);
    title(TITLE);
    drawnow
end


