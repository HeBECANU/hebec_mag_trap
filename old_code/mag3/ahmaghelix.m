%%Calculates the b-field of a 1 helical turn antihelmholtz coil
%%configuration
function[Bmat]=ahmaghelix(k,R,z1,wirethickness,Axx,Ayy,Azz,intsteps); 

Bmat=zeros(intsteps,3);

thetastart=2*pi;thetastop=0;
zstart=z1/2;
[Rx,Ry,Rz,dtheta,theta]=makehelixvector(R,zstart,wirethickness,thetastart,thetastop,intsteps); 
Bmat=helixturn(k,Rx,Ry,Rz,zstart,wirethickness,dtheta,theta,Axx,Ayy,Azz);

thetastart=0;thetastop=2*pi;
zstart=-z1/2;
[Rx,Ry,Rz,dtheta,theta]=makehelixvector(R,zstart,wirethickness,thetastart,thetastop,intsteps); 

Bmat=Bmat+helixturn(k,Rx,Ry,Rz,zstart,wirethickness,dtheta,theta,Axx,Ayy,Azz);
Bmat=sum(Bmat,1);





















