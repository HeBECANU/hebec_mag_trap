%%%%%%% subroutine for calculating magnetic field components at (Ax,Ay,Az)
% for one 'leaf' of a cloverleaf configuration with one loop.

function[Bmat]=coilmag(k,R,z1,Axx,Ayy,Azz,intsteps); 

Bmat=zeros(intsteps,3);

z=z1+1e-12;

thetastart=2*pi;thetastop=0;


[Rx,Ry,Rz,dtheta,theta]=makebiasvector(R,z,thetastart,thetastop,intsteps); 

  Bmat=bias(k,Rx,Ry,Rz,dtheta,theta,Axx,Ayy,Azz);
  
 
  Bmat=sum(Bmat,1);











