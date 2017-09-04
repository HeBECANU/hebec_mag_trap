%calculates the magnetic field for less than or equal to one helical turn
% offset in the x-dir

function[Btotal]=helixturnoff(k,offset,Rx,Ry,Rz,zstart,wirethickness,dtheta,theta,Axx,Ayy,Azz) 

% we assume that the wirethickness dictates the pitch and that the wires
% are snug up against each other.

omega=pi/wirethickness; %angular rotation frequency of the helix
b=omega/(2*pi)*wirethickness; % z-motion of the helix

Btotal=magarcmathelixoff(k,offset,Rx,Ry,Rz,zstart,Axx,Ayy,Azz,dtheta,theta,b,omega);


