%%%%%%% subroutine for calculating magnetic field components at (Ax,Ay,Az)

function[Bmat]=curvmagextra(k,R,z1,Axx,Ayy,Azz,intsteps,phi); 

Bmat=zeros(intsteps,3);

z=-z1/2+1e-12;thetastart=phi;thetastop=0; 
[Rx,Ry,Rz,dtheta,theta]=makebiasvector(R,z,thetastart,thetastop,intsteps); 

  Bmat=bias(k,Rx,Ry,Rz,dtheta,theta,Axx,Ayy,Azz);
  
  z=z1/2+1e-12;thetastart=pi+phi;thetastop=pi; 
  [Rx,Ry,Rz,dtheta,theta]=makebiasvector(R,z,thetastart,thetastop,intsteps); 

  Bmat=Bmat+bias(k,Rx,Ry,Rz,dtheta,theta,Axx,Ayy,Azz);

  Bmat=sum(Bmat,1);











