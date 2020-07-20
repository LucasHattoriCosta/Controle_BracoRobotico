G = tf(1,[1 1 0]);
Gc = 5;
T = [0:0.2:10];
[y,t] = step(feedback(G*Gc,1),T); 
[u,t] = step(feedback(Gc,G),T); 

plot(t,y,t,u);
xlabel('Tempo');
legend('y','u')

f=figure;
ax1 = axes('Parent',f);
l = line([0 5],[0 0]);
ht = text(5,6,'t=0');
set(l,'Linewidth',4);
axis equal;
set(ax1,'XLim',[-1 7]);
set(ax1,'YLim',[-1 7]);

ax2 = axes('Parent',f,'Units','normalized','Position',[0.75 0.2 0.1 0.3]);
b = bar(1)
set(ax2,'XLim',[-0 2])
set(ax2,'YLim',[-5 10])

set(f,'CurrentAxes',ax1);

for i=1:4
    ang = 30 * (i-1);
    l_grid(i) = line([0 10*cosd(ang)],[0 10*sind(ang)]);
    set(l_grid(i),'LineStyle','--','Color','k');
end;

v = VideoWriter('newfile.avi');
open(v)
for i=1:length(t)
    ang = 60 * y(i);
    set(l,'XData',[0 5*cosd(ang)]);
    set(l,'YData',[0 5*sind(ang)]);
    set(ht,'String',['t=',num2str(t(i),2)],'Visible','on');
    set(b,'YData',u(i)*5);
    frame = getframe(gcf);
    writeVideo(v,frame)
end;
close(v)

    

 