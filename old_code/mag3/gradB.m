%%%function that finds the GRAD of a B-field in G/m
%% in this case antihelquic gives the b-field
function[gb] = gradB(x,y,z,step)

x1=x-step;
x2=x+step;
y1=y-step;
y2=y+step;
z1=z-step;
z2=z+step;

bx1=aqs(x1,y,z);
bx2=aqs(x2,y,z);
dbx=(bx1(1)-bx2(1))/(2*step);

by1=aqs(x,y1,z);
by2=aqs(x,y2,z);
dby=(by1(1)-by2(1))/(2*step);

bz1=aqs(x,y,z1);
bz2=aqs(x,y,z2);
dbz=(bz1(1)-bz2(1))/(2*step);

gb(1)=dbx;
gb(2)=dby;
gb(3)=dbz;

gb=gb/1e-4;