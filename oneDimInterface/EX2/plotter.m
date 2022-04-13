clear, clf
figure(1)
S1 = dir('solution1*.ext');
S2 = dir('solution2*.ext');
N = length(S1);
s1 = load(S1(1).name);
n1 = length(s1) - 1;
h1 = 2.0/n1;
x1 = -1.0 + (0:n1)'*h1;

s2 = load(S2(1).name);
n2 = length(s2) - 1;
h2 = 2.0/n2;
x2 = -1.0 + (0:n2)'*h2;
% vidObj = VideoWriter('interface.mp4','MPEG-4');
% vidObj.FrameRate = 32;
% open(vidObj);

C = 1;
D = 2;
T = 2*C/(D+C);
R = T-1;
for i = 1:N
    s1 = load(S1(i).name);
    s2 = load(S2(i).name);
    plot(x1,s1./(x1<=0.0),x2,s2./(x2>=0.0),'LineWidth',2)
    axis([-1 1 -1 1])
    drawnow
    %     currFrame = getframe(gcf);
    %     writeVideo(vidObj,currFrame);
end
% close(vidObj)
Tend = 0.8
u2ex = T*0.5*exp(-500*((x2-D*0.3)/D).^2)./(x2>=0.0);
u1ex = R*0.5*exp(-500*((x1+C*0.3)/C).^2)./(x1<=0.0) - 0.5*exp(-500*((x1+C*0.7)/C).^2)./(x1<=0.0);
plot(x1,s1./(x1<=0.0),x2,s2./(x2>=0.0),...
     x2,u2ex,'--',...
     x1,u1ex,'--','LineWidth',2)

