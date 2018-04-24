%%%%%%% subroutine for calculating magnetic field components at (Ax,Ay,Az)

function[Bmat]=Qmagy(k,R,y1,Axx,Ayy,Azz,intsteps); 

Bmat=zeros(intsteps,3);

y=-y1/2+1e-12;thetastart=0;thetastop=2*pi; 

[Rx,Ry,Rz,dtheta,theta]=makeQvectory(R,y,thetastart,thetastop,intsteps); 

  Bmat=biasy(k,Rx,Ry,Rz,dtheta,theta,Axx,Ayy,Azz);
  
  Bmat=sum(Bmat,1);
  











