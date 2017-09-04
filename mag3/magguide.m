clear all
close all

x1=-1/1000;
y1=1/1000;
x2=1/1000;
y2=1/1000;
x3=-1/1000;
y3=-1/1000;
x4=1/1000;
y4=-1/1000;

mu0=(4*pi*1e-7);
i=1;
const=mu0*i/(2*pi)/1e-4;

range=.1/1000;
xsteps=100;
ysteps=100;
ys=0;
xs=0;
y=ys-range:range/ysteps:ys+range;
x=xs-range:range/xsteps:xs+range;

% first along x

r1i=1./sqrt((x1-x).^2+(y1-ys).^2);
r2i=1./sqrt((x2-x).^2+(y2-ys).^2);
r3i=1./sqrt((x3-x).^2+(y3-ys).^2);
r4i=1./sqrt((x4-x).^2+(y4-ys).^2);

B=const*(r1i-r2i);%-r3i+r4i);
B=sqrt(B.^2);
plot(x,B)
