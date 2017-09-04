% program that calculates the B-field of a 4-wire guide
% the wires are located at (x,y)= (XX,YY)
function[B]=fourwirehumpfield(x)

Mu0 = 4*pi*1e-7;
I=10;
Bc=Mu0*I/(2*pi)/1e-4; %put in Gauss

XX=.30/100;
YY=.30/100;
y=0;
yp=y+YY;
ym=y-YY;
xp=x+XX;
xm=x-XX;

Bw1ppx=yp./(xp.^2+yp.^2);
Bw1ppy=xp./(xp.^2+yp.^2);
Bw2pmx=ym./(xp.^2+ym.^2);
Bw2pmy=xp./(xp.^2+ym.^2);
Bw3mmx=ym./(xm.^2+ym.^2);
Bw3mmy=xm./(xm.^2+ym.^2);
Bw4mpx=yp./(xm.^2+yp.^2);
Bw4mpy=xm./(xm.^2+yp.^2);
    
    BBx= Bw1ppx+Bw2pmx-Bw3mmx-Bw4mpx;
     BBy= Bw1ppy+Bw2pmy-Bw3mmy-Bw4mpy;
    
 B=Bc*sqrt(BBx.^2+BBy.^2);
