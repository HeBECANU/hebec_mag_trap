%%%%%%% subroutine for calculating magnetic field components at (Ax,Ay,Az)
% for one 'leaf' of a cloverleaf configuration with one loop.

function[Bmat]=extramag(k,Rin,Rout,z1,Axx,Ayy,Azz,intsteps); 

Bmat=zeros(intsteps,3);

% extra lengths for one cloverleaf coil
 z=-z1/2+1e-12;thetainstart=pi;,thetainstop=pi/2;,thetaoutstart=pi/2;,thetaoutstop=0;
 [Rxi,Ryi,Rzi,Rxo,Ryo,Rzo,dtheta1,dtheta2,theta1,theta2]=makeextravector(Rin,Rout,...
    z,thetaoutstart,thetaoutstop,thetainstart,thetainstop,intsteps); 

  Bmat=extras(k,Rxi,Ryi,Rzi,Rxo,Ryo,Rzo,dtheta1,dtheta2,theta1,...
     theta2,Axx,Ayy,Azz);
  
  z=-z1/2+1e-12;thetainstart=2*pi;,thetainstop=3*pi/2;,thetaoutstart=3*pi/2;,thetaoutstop=pi;
 [Rxi,Ryi,Rzi,Rxo,Ryo,Rzo,dtheta1,dtheta2,theta1,theta2]=makeextravector(Rin,Rout,...
    z,thetaoutstart,thetaoutstop,thetainstart,thetainstop,intsteps); 

  Bmat=Bmat+extras(k,Rxi,Ryi,Rzi,Rxo,Ryo,Rzo,dtheta1,dtheta2,theta1,...
     theta2,Axx,Ayy,Azz);
  
  % other cloverleaf coil
  
   z=z1/2+1e-12;thetainstart=3*pi/2;,thetainstop=2*pi;,thetaoutstart=pi;,thetaoutstop=3*pi/2;
 [Rxi,Ryi,Rzi,Rxo,Ryo,Rzo,dtheta1,dtheta2,theta1,theta2]=makeextravector(Rin,Rout,...
    z,thetaoutstart,thetaoutstop,thetainstart,thetainstop,intsteps); 

  Bmat=Bmat+extras(k,Rxi,Ryi,Rzi,Rxo,Ryo,Rzo,dtheta1,dtheta2,theta1,...
     theta2,Axx,Ayy,Azz);
  
  z=z1/2+1e-12;thetainstart=pi/2;,thetainstop=pi;,thetaoutstart=0;,thetaoutstop=pi/2;
 [Rxi,Ryi,Rzi,Rxo,Ryo,Rzo,dtheta1,dtheta2,theta1,theta2]=makeextravector(Rin,Rout,...
    z,thetaoutstart,thetaoutstop,thetainstart,thetainstop,intsteps); 

  Bmat=Bmat+extras(k,Rxi,Ryi,Rzi,Rxo,Ryo,Rzo,dtheta1,dtheta2,theta1,...
     theta2,Axx,Ayy,Azz);


  Bmat=sum(Bmat,1);












