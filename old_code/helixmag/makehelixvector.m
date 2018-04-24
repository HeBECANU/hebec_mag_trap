%%%%%%% subroutine for making the relevant coordinates for a <=1 turn helix
function[Rx,Ry,Rz,dtheta,theta,b,omega]=makehelixvector(R,zstart,wirethickness,thetastart,thetastop,intsteps) 

% now arc coordinates
% remember to get the coordinates right in terms of direction of the current
% we assume that the wirethickness dictates the pitch and that the wires
% are snug up against each other.

omega=pi/wirethickness; %angular rotation frequency of the helix
b=omega/(2*pi)*wirethickness; % z-motion of the helix

dtheta=(thetastop-thetastart)/intsteps;
theta = thetastart+dtheta/2:dtheta:thetastop-dtheta/2;


helixdtheta=abs(thetastop-thetastart)/intsteps;
helixtheta=helixdtheta/2:helixdtheta:abs(thetastop-thetastart)-helixdtheta/2;

    
    
t=helixtheta/omega; %parametric parameter for the helix we use a seperate theta here to ensure
                    % the helix spirals in the same way independent of the
                    % given theta
   
   Rx=R*cos(theta);
   Ry=R*sin(theta);
   if zstart>0
   Rz=zstart+(b*t);
   end
   if zstart<0
   Rz=zstart-(b*t);
   end
