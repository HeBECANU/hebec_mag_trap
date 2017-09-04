%%%%%%% subroutine for making the relevant cloverleaf coordinates
function[Rxi,Ryi,Rzi,Rxo,Ryo,Rzo,dtheta1,dtheta2,theta1,theta2]=...
   makeextravector(Rin,Rout,z1stop,thetaoutstart,thetaoutstop,thetainstart,thetainstop...
,intsteps) 

% now arc coordinates
% remember to get the coordinates right in terms of direction of the current

dtheta1=(thetainstop-thetainstart)/intsteps;
theta1 = thetainstart+dtheta1/2:dtheta1:thetainstop-dtheta1/2;
   
   Rxi=Rin*cos(theta1);
   Ryi=Rin*sin(theta1);
   Rzi=z1stop;
   
dtheta2=(thetaoutstop-thetaoutstart)/intsteps;
theta2 = thetaoutstart+dtheta2/2:dtheta2:thetaoutstop-dtheta2/2;
  
   Rxo=Rout*cos(theta2);
   Ryo=Rout*sin(theta2);
   Rzo=z1stop;

