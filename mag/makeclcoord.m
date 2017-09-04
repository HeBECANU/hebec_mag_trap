%%%%%%% subroutine for making the relevant cloverleaf coordinates
function[startx1,starty1,startz1,stopx1,stopy1,...
   stopz1,startx2,starty2,startz2,stopx2,stopy2,stopz2,Rxi,Ryi,...
   Rzi,Rxo,Ryo,Rzo,dtheta,thetain,thetaout]=...
makevector(x1start,x1stop,y1start,y1stop,z1start,z1stop,x2start,x2stop,y2start,...
y2stop,z2start,z2stop,Rin,Rout,thetaoutstart,thetaoutstop,thetainstart,thetainstop...
,intsteps) 

%first define clover leaf coordinates
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

% now arc coordinates
% remember to get the coordinates right in terms of direction of the current

dtheta=(thetainstop-thetainstart)/intsteps;
theta1 = thetainstart+dtheta/2:dtheta/2:thetainstop-detheta/2;
   
   Rxi=Rin*cos(theta1);
   Ryi=Rin*sin(theta1);
   Rzi=z1stop;
   
dtheta=(thetaoutstop-thetaoutstart)/intsteps;
theta2 = thetaoutstart+dtheta/2:dtheta/2:thetaoutstop-detheta/2;
  
   Rxo=Rout*cos(theta2);
   Ryo=Rout*sin(theta2);
   Rzo=z1stop;

