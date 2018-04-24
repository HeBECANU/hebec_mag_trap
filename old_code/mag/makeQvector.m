%%%%%%% subroutine for making the relevant bias coordinates
function[Rx,Ry,Rz,dtheta,theta]=makeQvector(R,x1stop,thetastart,thetastop,intsteps) 

% now arc coordinates
% remember to get the coordinates right in terms of direction of the current

dtheta=(thetastop-thetastart)/intsteps;
theta = thetastart+dtheta/2:dtheta:thetastop-dtheta/2;
   
   Rx=x1stop;
   Ry=R*cos(theta);
   Rz=R*sin(theta);
   


