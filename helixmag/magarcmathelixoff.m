%%%%%%% subroutine for calculating magnetic field components of an arc at a
% particular distance.
% same as magarcmat accept it allows for motion in the z-direction
% i.e a helix
% offset in the x-direction
function[B]=magarcmathelixoff(k,offset,Rx,Ry,Rz,zstart,Ax,Ay,Az,dtheta,theta,b,omega) 

% ok let the origin be in the center of the trapping region

% the relevant distances are then:

rx=Rx-Ax+offset; %We are off axis in the x-axis
ry=Ry-Ay;
rz=Rz-Az;

% thus the vector r is given by

r=zeros(length(theta),3);
R=zeros(length(theta),3);
dlunitXrunit=zeros(length(theta),3);

r(:,1)=rx';
r(:,2)=ry';
r(:,3)=rz';

% and the magnitude of r is

rmag =sqrt(sum(r.^2,2));
rmag(:,2)=rmag(:,1);
rmag(:,3)=rmag(:,1);

% the unit vector in the r direction
runit=r./rmag;

% work out radius
R(:,1)= Rx';
R(:,2)= Ry';

Rmag =sqrt(sum(R.^2,2));
Rmag(:,2)=Rmag(:,1);
Rmag(:,3)=Rmag(:,1);

% now dl = dx (i) + dy (j) +dz (k)

dx=-1*Rmag(1,1).*sin(theta')*dtheta;
dy=Rmag(1,1).*cos(theta')*dtheta;
dz=b/omega*dtheta;

dl(:,1)=dx;
dl(:,2)=dy;
dl(:,3)=dz;

dlmag=sqrt(sum(dl.^2,2));
dlmag(:,2)=dlmag(:,1);
dlmag(:,3)=dlmag(:,1);

dlunit=dl./dlmag;

dlunitXrunit(:,1)=(dlunit(:,2).*runit(:,3))-(dlunit(:,3).*runit(:,2));
dlunitXrunit(:,2)=(dlunit(:,3).*runit(:,1))-(dlunit(:,1).*runit(:,3));
dlunitXrunit(:,3)=(dlunit(:,1).*runit(:,2))-(dlunit(:,2).*runit(:,1));

B=k*dlmag.*dlunitXrunit./rmag.^2;
