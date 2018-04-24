%%%%%%% subroutine for calculating magnetic field components at (Ax,Ay,Az)
% for one 'leaf' of a cloverleaf configuration with one loop.
function[Bmat]=ahmagextra(k,R,thetastart,thetastop,z1,Axx,Ayy,Azz,intsteps); 

Bmat=zeros(intsteps,3);
z=z1/2+1e-12;

[Rx,Ry,Rz,dtheta,theta]=makebiasvector(R,z,thetastart,thetastop,intsteps); 
Bmat=bias(k,Rx,Ry,Rz,dtheta,theta,Axx,Ayy,Azz);
z=-z1/2+1e-12;
thetastart2=thetastop;thetastop2=thetastart;

[Rx,Ry,Rz,dtheta,theta]=makebiasvector(R,z,thetastart2,thetastop2,intsteps); 

Bmat=Bmat+bias(k,Rx,Ry,Rz,dtheta,theta,Axx,Ayy,Azz);
Bmat=sum(Bmat,1);





















