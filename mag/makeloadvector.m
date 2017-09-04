%%%%%%% subroutine for making the relevant loading guide coordinates
function[startx1,starty1,startz1,stopx1,stopy1,stopz1,...
        startx2,starty2,startz2,stopx2,stopy2,stopz2,startx3,starty3,startz3,stopx3,stopy3,stopz3,...
        startx4,starty4,startz4,stopx4,stopy4,stopz4]=makeloadvector(x1start,y1start,z1start,x1stop,y1stop,z1stop,...
        x2start,y2start,z2start,x2stop,y2stop,z2stop,x3start,y3start,z3start,x3stop,y3stop,z3stop,...
        x4start,y4start,z4start,x4stop,y4stop,z4stop,intsteps) 

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

dx3=(x3stop-x3start)/intsteps;
dy3=(y3stop-y3start)/intsteps;
dz3=(z3stop-z3start)/intsteps;
startx3=x3start:dx3:x3stop-dx3;
starty3=y3start:dy3:y3stop-dy3;
startz3=z3start:dz3:z3stop-dz3;
stopx3=x3start+dx3:dx3:x3stop;
stopy3=y3start+dy3:dy3:y3stop;
stopz3=z3start+dz3:dz3:z3stop;

dx4=(x4stop-x4start)/intsteps;
dy4=(y4stop-y4start)/intsteps;
dz4=(z4stop-z4start)/intsteps;
startx4=x4start:dx4:x4stop-dx4;
starty4=y4start:dy4:y4stop-dy4;
startz4=z4start:dz4:z4stop-dz4;
stopx4=x4start+dx4:dx4:x4stop;
stopy4=y4start+dy4:dy4:y4stop;
stopz4=z4start+dz4:dz4:z4stop;


