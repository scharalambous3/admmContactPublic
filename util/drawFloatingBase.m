function drawFloatingBase(hax, c, f, xDes, rc1, rc2)
W = 0.75;  % cart width
H = .5; % cart height
figure(3)
cla
%plot(hax, [-(c(1)+W/2) (c(1)+W/2)],[-(c(2)+H/2) (c(2)+H/2)],'k','LineWidth',2)
hold on
%for i = 1:length(t)
%x = xTraj(1,i);
% plot(hax, [-(xDes+W/2) (xDes+W/2)],[0 0],'k','LineWidth',2)
%hold on
%rectangle('Position',[c(1)-W/2,c(2)-H/2,W,H],'Curvature',.1)

x = [c(1) - W/2, c(1) + W/2];
y = [c(2) - H/2, c(2) + H/2];
% h = plot([x(1) x(2) x(2) x(1) x(1)],...
%    [y(1) y(1) y(2) y(2) y(1)]);
%polyin = fill(X,Y,'k');
%rotate(h,[0 0 1], rad2deg(theta), c)
theta = 0;
R = [cos(theta) -sin(theta); sin(theta) cos(theta)];
C  = repmat(c, 1, 5);

V = [x(1), x(2), x(2), x(1), x(1); y(1), y(1), y(2), y(2), y(1)];
%V = get(polyin,'Vertices')';  % get the current set of vertices
V = R*(V - C) + C;             % do the rotation relative to the centre of the square
%set(polyin,'Vertices',V'); 
x = V(1,:);
y = V(2,:);

plot(x, y);
yline(0.5, '--');
yline(0.0, '-');

quiver(rc1(1), rc1(2), f(1)/100, f(2)/100)
quiver(rc2(1), rc2(2), f(3)/100, f(4)/100)

paddingX=1.25;
% set(gca,'YTick',[])
% set(gca,'XTick',[])
xlim([-W xDes+(W/2)*paddingX]);
ylim([-5 15]);
set(gcf,'Position',[100 550 5000 5000])
% box off

%dt_viz=0.1
drawnow
%gif
pause
%pause(dt_viz)
axis equal
    
%end
hold off