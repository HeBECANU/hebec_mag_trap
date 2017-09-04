%%%%%%% subroutine for calculating magnetic field components of an arc at a
% particular distance.
function[B]=magarcmatoffaxis(k,offset,Rx,Ry,Rz,Ax,Ay,Az,dtheta,theta) 

% the relevant distances are then:
rx=Rx-Ax;% +offset if you want to be off in both x and y
ry=Ry-Ay+offset;
rz=Rz-Az;

% thus the vector r is given by
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
R(:,1)= Rx';
R(:,2)=Ry';
Rmag =sqrt(sum(R.^2,2));
Rmag(:,2)=Rmag(:,1);
Rmag(:,3)=Rmag(:,1);

% now have to work out the direction of dl
dy=Rmag(1,1).*cos(theta')*dtheta;
dx=-1*Rmag(1,1).*sin(theta')*dtheta;
dl(:,1)=dx;
dl(:,2)=dy;
dl(:,3)=0;
dlmag=sqrt(sum(dl.^2,2));
dlmag(:,2)=dlmag(:,1);
dlmag(:,3)=dlmag(:,1);

dlunit=dl./dlmag;
dlmag=abs(Rmag*dtheta);
dlunitXrunit=zeros(length(Rx),3);
dlunitXrunit(:,1)=(dlunit(:,2).*runit(:,3))-(dlunit(:,3).*runit(:,2));
dlunitXrunit(:,2)=(dlunit(:,3).*runit(:,1))-(dlunit(:,1).*runit(:,3));
dlunitXrunit(:,3)=(dlunit(:,1).*runit(:,2))-(dlunit(:,2).*runit(:,1));
B=k*dlmag.*dlunitXrunit./rmag.^2;