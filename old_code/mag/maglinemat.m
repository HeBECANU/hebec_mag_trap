%%%%%%% subroutine for calculating magnetic field components at (Ax,Ay,Az)
% of a line which starts at (startx,starty,startz)
% and extends dl distance.

function[B]=maglinemat(k,startx,starty,startz,stopx,stopy,stopz,Ax,Ay,Az) 

% ok let the origin be in the center of the trapping region

% the relevant distances are then:
rx=(startx+(stopx-startx)/2)-Ax;
ry=(starty+(stopy-starty)/2)-Ay;
rz=(startz+(stopz-startz)/2)-Az;

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

% now calculate dl 

% now have to work out the direction of dl
dy=(stopy-starty);
dx=(stopx-startx);
dz=(stopz-startz);
dl(:,1)=dx';
dl(:,2)=dy';
dl(:,3)=dz';
dlmag=sqrt(sum(dl.^2,2));
dlmag(:,2)=dlmag(:,1);
dlmag(:,3)=dlmag(:,1);

dlunit=dl./dlmag;


dlunitXrunit=zeros(length(startx),3);
dlunitXrunit(:,1)=(dlunit(:,2).*runit(:,3))-(dlunit(:,3).*runit(:,2));
dlunitXrunit(:,2)=(dlunit(:,3).*runit(:,1))-(dlunit(:,1).*runit(:,3));
dlunitXrunit(:,3)=(dlunit(:,1).*runit(:,2))-(dlunit(:,2).*runit(:,1));


B=k*dlmag.*dlunitXrunit./rmag.^2;
