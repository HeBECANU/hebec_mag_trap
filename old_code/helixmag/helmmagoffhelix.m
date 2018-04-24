%%%%%%% subroutine for calculating magnetic field components at (Ax,Ay,Az)
% calculates the field due to a anti-helmholtz helical coil offset in the
% x-axis
function[Bmat]=helmmagoffhelix(k,offset,R,thetastart,thetastop,z1,wirethickness,Axx,Ayy,Azz,intsteps); 

Bmat=zeros(intsteps,3);
zstart=z1/2;
[Rx,Ry,Rz,dtheta,theta]=makehelixvector(R,zstart,wirethickness,thetastart,thetastop,intsteps); 
Bmat=helixturnoff(k,offset,Rx,Ry,Rz,zstart,wirethickness,dtheta,theta,Axx,Ayy,Azz);

zstart=-z1/2;

[Rx,Ry,Rz,dtheta,theta]=makehelixvector(R,zstart,wirethickness,thetastart,thetastop,intsteps); 
Bmat=Bmat+helixturnoff(k,offset,Rx,Ry,Rz,zstart,wirethickness,dtheta,theta,Axx,Ayy,Azz);
Bmat=sum(Bmat,1);

