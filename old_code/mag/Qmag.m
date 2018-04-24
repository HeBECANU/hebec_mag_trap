%%%%%%% subroutine for calculating magnetic field components at (Ax,Ay,Az)

function[Bmat]=Qmag(k,R,x1,Axx,Ayy,Azz,intsteps); 

Bmat=zeros(intsteps,3);

x=x1/2+1e-12;thetastart=0;thetastop=2*pi; 

[Rx,Ry,Rz,dtheta,theta]=makeQvector(R,x,thetastart,thetastop,intsteps); 

  Bmat=biasx(k,Rx,Ry,Rz,dtheta,theta,Axx,Ayy,Azz);
  
  Bmat=sum(Bmat,1);
  











