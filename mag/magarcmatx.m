%%%%%%% subroutine for calculating magnetic field components of an arc at a
% particular distance.

function[B]=magarcmatx(k,Rx,Ry,Rz,Ax,Ay,Az,dtheta,theta) 

% ok let the origin be in the center of the trapping region

% the relevant distances are then:

rx=Rx-Ax;
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

% now dl = R * dtheta
% Where R=magnitude of the radius

R(:,1)= Rz';
R(:,2)= Ry';
Rmag =sqrt(sum(R.^2,2));
Rmag(:,2)=Rmag(:,1);
Rmag(:,3)=Rmag(:,1);

% now have to work out the direction of dl

dz=Rmag(1,1).*cos(theta')*dtheta;
dy=-1*Rmag(1,1).*sin(theta')*dtheta;

dl(:,2)=dy;
dl(:,3)=dz;
dl(:,1)=0;

dlmag=sqrt(sum(dl.^2,2));
dlmag(:,2)=dlmag(:,1);
dlmag(:,3)=dlmag(:,1);

dlunit=dl./dlmag;
dlmag=abs(Rmag*dtheta);

dlunitXrunit(:,1)=(dlunit(:,2).*runit(:,3))-(dlunit(:,3).*runit(:,2));
dlunitXrunit(:,2)=(dlunit(:,3).*runit(:,1))-(dlunit(:,1).*runit(:,3));
dlunitXrunit(:,3)=(dlunit(:,1).*runit(:,2))-(dlunit(:,2).*runit(:,1));

B=k*dlmag.*dlunitXrunit./rmag.^2;
