u = dir('u*.ext');
S = size(load(u(1).name));
k31 = 6.3801618959239383506237;
k33 = 13.01520072169843441983;
gamma = k31/k33;
alpha = 1-gamma;
%% Set reference element
rs = 0.0;
re = 1.0;
ss = 0.0;
se = 1.0;
ns = S(1)-1;
nr = S(2)-1;
hr = (re-rs)/nr;
hs = (se-ss)/ns;
r = rs + (0:nr)*hr;
s = ss + (0:ns)*hs;

%% Map to polar

x = zeros(nr+1,ns+1);
for i = 0:nr
    for j = 0:ns
        x(i+1,j+1) = (alpha*r(1+i)+gamma)*cos(2*pi*s(j+1));
    end
end

y = zeros(nr+1,ns+1);
for i = 0:nr
    for j = 0:ns
        y(i+1,j+1) = (alpha*r(1+i)+gamma)*sin(2*pi*s(j+1));
    end
end
 
% v = VideoWriter('polarWave.mp4','MPEG-4');
% v.FrameRate = 10;
% open(v);
for n = 1:length(u)
    uplot = load(u(n).name);
    surf(x,y,uplot');
    shading interp
    
    drawnow
     % Capture the plot as an image 
%      frame = getframe(gcf); 
%      writeVideo(v,frame);
end 
% close(v);

      