clear all
y=1:8
x=(1:8)*.3
plot(x,y)
H=gca
maxx=str2num(get(H,'xticklabel'));
scalex=(maxx(length(maxx))-maxx(length(maxx)-1))
maxx=max(maxx)
maxy=str2num(get(H,'yticklabel'));
scaley=(maxy(length(maxy))-maxy(length(maxy)-1))
maxy=max(maxy)