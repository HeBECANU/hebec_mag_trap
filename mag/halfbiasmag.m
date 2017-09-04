%%%%%%% subroutine for calculating magnetic field components at (Ax,Ay,Az)
%%%  calculates the extra due to extra half wrap
% note that the half wraps are not symmetric

function[Bmat]=halfbiasmag(k,R,z1,Axx,Ayy,Azz,intsteps); 

Bmat=zeros(intsteps,3);

z=-z1/2+1e-12;thetastart=0;thetastop=pi; 
[Rx,Ry,Rz,dtheta,theta]=makebiasvector(R,z,thetastart,thetastop,intsteps); 

  Bmat=bias(k,Rx,Ry,Rz,dtheta,theta,Axx,Ayy,Azz);
  
  z=z1/2+1e-12;thetastart=pi;thetastop=2*pi; 
[Rx,Ry,Rz,dtheta,theta]=makebiasvector(R,z,thetastart,thetastop,intsteps); 

  Bmat=Bmat+bias(k,Rx,Ry,Rz,dtheta,theta,Axx,Ayy,Azz);

  Bmat=sum(Bmat,1);











