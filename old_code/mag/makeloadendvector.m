%%%%%%% subroutine for making the relevant loading guide coordinates
function[startx1,starty1,startz1,stopx1,stopy1,stopz1,...
        startx2,starty2,startz2,stopx2,stopy2,stopz2]=makeloadendvector(x1start,y1start,z1start,x1stop,y1stop,z1stop,...
        x2start,y2start,z2start,x2stop,y2stop,z2stop,intsteps) 

% define line coordinates first
dx1=(x1stop-x1start)/intsteps;
dy1=(y1stop-y1start)/intsteps;
dz1=(z1stop-z1start)/intsteps;
startx1=x1start:dx1:x1stop-dx1;
starty1=y1start:dy1:y1stop-dy1;
startz1=z1start:dz1:z1stop-dz1;
stopx1=x1start+dx1:dx1:x1stop;
stopy1=y1start+dy1:dy1:y1stop;
stopz1=z1start+dz1:dz1:z1stop;

dx2=(x2stop-x2start)/intsteps;
dy2=(y2stop-y2start)/intsteps;
dz2=(z2stop-z2start)/intsteps;
startx2=x2start:dx2:x2stop-dx2;
starty2=y2start:dy2:y2stop-dy2;
startz2=z2start:dz2:z2stop-dz2;
stopx2=x2start+dx2:dx2:x2stop;
stopy2=y2start+dy2:dy2:y2stop;
stopz2=z2start+dz2:dz2:z2stop;


