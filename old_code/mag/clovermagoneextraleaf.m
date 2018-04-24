%%%%%%% subroutine for calculating magnetic field components at (Ax,Ay,Az)
% for one 'leaf' of a cloverleaf configuration with one loop.

function[Bmat]=clovermagoneextraleaf(k,Rin,Rout,z1,Axx,Ayy,Azz,intsteps,symetry,extra); 

   Bmat=zeros(intsteps,3);

   % cloverleaf 1
  x1start=-Rout;,x1stop=-Rin;,y1start=0;,y1stop=1e-20;,z1start=-z1/2;,z1stop=-z1/2+1e-12;
  x2start=0;,x2stop=1e-20;,y2start=Rin;,y2stop=Rout;,z2start=-z1/2;,z2stop=-z1/2+1e-12;
  thetainstart=pi;,thetainstop=pi/2;,thetaoutstart=pi/2;,thetaoutstop=pi;
[startx1,starty1,startz1,stopx1,stopy1,stopz1,startx2,starty2,startz2,stopx2,stopy2,...
      stopz2,Rxi,Ryi,Rzi,Rxo,Ryo,Rzo,dtheta1,dtheta2,theta1,theta2]=makevector(x1start,...
   x1stop,y1start,y1stop,z1start,z1stop,x2start,x2stop,y2start,y2stop,z2start,z2stop,...
Rin,Rout,thetaoutstart,thetaoutstop,thetainstart,thetainstop,intsteps); 

  Bmat=clmagmat(k,startx1,starty1,startz1,stopx1,stopy1,...
   stopz1,startx2,starty2,startz2,stopx2,stopy2,stopz2,Rxi,Ryi,...
   Rzi,Rxo,Ryo,Rzo,dtheta1,dtheta2,theta1,theta2,Axx,Ayy,Azz,extra);

Bmat=sum(Bmat,1);












