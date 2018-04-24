%%%%%%% subroutine for calculating magnetic field components at (Ax,Ay,Az)
% for one 'leaf' of a cloverleaf configuration with one loop.

function[Btotal]=clmagmat(k,startx1,starty1,startz1,stopx1,stopy1,...
   stopz1,startx2,starty2,startz2,stopx2,stopy2,stopz2,Rxi,Ryi,...
   Rzi,Rxo,Ryo,Rzo,dthetain,dthetaout,thetain,thetaout,Ax,Ay,Az,extra) 

% first calculate the magfield due to the line elements
% we assumme the user has input the start and stop values so that
% the direction of the current flows from the start to the stop coordinate.

Bline1=maglinemat(k,startx1,starty1,startz1,stopx1,stopy1,stopz1,Ax,Ay,Az);
Bline2=maglinemat(k,startx2,starty2,startz2,stopx2,stopy2,stopz2,Ax,Ay,Az);

%ok now for the arcs once again we assume the user entered the 
% right direction

Binarc=magarcmat(k,Rxi,Ryi,Rzi,Ax,Ay,Az,dthetain,thetain);

if extra <1 
Boutarc=magarcmat(k,Rxo,Ryo,Rzo,Ax,Ay,Az,dthetaout,thetaout);
else Boutarc=0;
end

% now lets add all the fields up
Btotal=Bline1+Bline2+Binarc+Boutarc;

