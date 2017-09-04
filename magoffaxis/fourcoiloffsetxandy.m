%%%%%%% subroutine for calculating magnetic field components at (Ax,Ay,Az)
% for one 'leaf' of a cloverleaf configuration with one loop.
function[Bmat]=fourcoiloffsetxandy(k,offsetx,offsety,R,z1,Axx,Ayy,Azz,intsteps); 

Bmat=zeros(intsteps,3);
z=z1/2+1e-12;
thetastart=2*pi;thetastop=0;
[Rx,Ry,Rz,dtheta,theta]=makebiasvector(R,z,thetastart,thetastop,intsteps); 
Bmat=biasoffxandy(k,offsetx,offsety,Rx,Ry,Rz,dtheta,theta,Axx,Ayy,Azz); %x=+ y=+

thetastart=0;thetastop=2*pi;
[Rx,Ry,Rz,dtheta,theta]=makebiasvector(R,z,thetastart,thetastop,intsteps); 
Bmat=Bmat+biasoffxandy(k,offsetx,-offsety,Rx,Ry,Rz,dtheta,theta,Axx,Ayy,Azz); %x=+ y=-

thetastart=2*pi;thetastop=0;
[Rx,Ry,Rz,dtheta,theta]=makebiasvector(R,z,thetastart,thetastop,intsteps); 
Bmat=Bmat+biasoffxandy(k,-offsetx,-offsety,Rx,Ry,Rz,dtheta,theta,Axx,Ayy,Azz); %x=- y=-

thetastart=0;thetastop=2*pi;
[Rx,Ry,Rz,dtheta,theta]=makebiasvector(R,z,thetastart,thetastop,intsteps); 
Bmat=Bmat+biasoffxandy(k,-offsetx,offsety,Rx,Ry,Rz,dtheta,theta,Axx,Ayy,Azz); %x=- y=+

Bmat=sum(Bmat,1);

