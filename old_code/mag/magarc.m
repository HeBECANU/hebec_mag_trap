%%%%%%% subroutine for calculating magnetic field components of an arc at a
% particular distance.

function[B]=magarc(k,Rx,Ry,Rz,Ax,Ay,Az,dtheta,theta) 

% ok let the origin be in the center of the trapping region

% the relevant distances are then:

rx=Rx-Ax;
ry=Ry-Ay;
rz=Rz-Az;

% thus the vector r is given by
r(1,1)=rx;
r(1,2)=ry;
r(1,3)=rz;

% and the magnitude of r is

rmag =sqrt(sum(r.^2));

% the unit vector in the r direction

runit=r./rmag;

% now dl = R * dtheta
% Where R=magnitude of the radius

R(1,1)= Rx;
R(1,2)=Ry;
Rmag =sqrt(sum(R.^2));

% now have to work out the direction of dl

dy=Rmag*cos(theta)*dtheta;
dx=-1*Rmag*sin(theta)*dtheta;

dl(1,1)=dx;
dl(1,2)=dy;
dl(1,3)=0;
dlmag=sqrt(sum(dl.^2));
dlunit=dl/dlmag;

dlmag=Rmag*dtheta;

B=k*dlmag*cross(dlunit,runit)/rmag^2;
