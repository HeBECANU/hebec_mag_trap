%%Calculates the b-field of a 1 helical turn antihelmholtz coil
%%configuration
function[Bmat]=ahmaghelixextra(k,R,thetastart,thetastop,z1,wirethickness,Axx,Ayy,Azz,intsteps); 

Bmat=zeros(intsteps,3);

zstart=z1/2;
[Rx,Ry,Rz,dtheta,theta]=makehelixvector(R,zstart,wirethickness,thetastart,thetastop,intsteps); 
Bmat=helixturn(k,Rx,Ry,Rz,zstart,wirethickness,dtheta,theta,Axx,Ayy,Azz);

thetastart2=thetastop;thetastop2=thetastart;
zstart=-z1/2;
[Rx,Ry,Rz,dtheta,theta]=makehelixvector(R,zstart,wirethickness,thetastart2,thetastop2,intsteps); 

Bmat=Bmat+helixturn(k,Rx,Ry,Rz,zstart,wirethickness,dtheta,theta,Axx,Ayy,Azz);
Bmat=sum(Bmat,1);





















