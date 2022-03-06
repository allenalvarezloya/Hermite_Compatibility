clear, clc
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
for i = 1:N
    s1 = load(S1(i).name);
    s2 = load(S2(i).name);
    plot(x1,s1,x2,s2,'LineWidth',2)
    axis([-1 1 -3 3])
    drawnow
%     currFrame = getframe(gcf);
%     writeVideo(vidObj,currFrame);
end
% close(vidObj)

