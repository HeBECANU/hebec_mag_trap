%%%%%%% subroutine for calculating magnetic field components at (Ax,Ay,Az)
% for one 'leaf' of a cloverleaf configuration with one loop.
function[Bmat]=onecoiloffsety(k,offset,R,z1,Axx,Ayy,Azz,intsteps); 

Bmat=zeros(intsteps,3);
z=z1/2+1e-12;
thetastart=2*pi;thetastop=0;
[Rx,Ry,Rz,dtheta,theta]=makebiasvector(R,z,thetastart,thetastop,intsteps); 
Bmat=biasoffy(k,offset,Rx,Ry,Rz,dtheta,theta,Axx,Ayy,Azz);

Bmat=sum(Bmat,1);

