function[]=key() 

H=gca;
xx=str2num(get(H,'xticklabel'));
minx=min(xx);
maxx=max(xx);
yy=str2num(get(H,'yticklabel'));
miny=min(yy);
maxy=max(yy);

xscaling=2;
yscaling=15;

xpos=maxx-((maxx-minx)/xscaling);
ypos=maxy-((maxy-miny)/yscaling);
l=text(xpos,ypos,'Axial');
set(l,'Fontsize',14);

gl=(get(l,'extent'));
gl=gl(3);
x=[(xpos+gl) maxx];
y=[ypos ypos];
line(x,y,'linestyle','-.','color','r')

   xpos=maxx-((maxx-minx)/xscaling);
ypos=maxy-(2*(maxy-miny)/yscaling);
l=text(xpos,ypos,'Radial');
get(l);
set(l,'Fontsize',14);

gl=(get(l,'extent'));
gl=gl(3);
x=[(xpos+gl) maxx];
y=[ypos ypos];
line(x,y,'linestyle','-','color','b')

xpos=maxx-((maxx-minx)/xscaling);
ypos=maxy-(3*(maxy-miny)/yscaling);
%l=text(xpos,ypos,'xy-axis');
gl=(get(l,'extent'));
gl=gl(3);
x=[(xpos+gl) maxx];
y=[ypos ypos];
%line(x,y,'linestyle',':','color','g')

