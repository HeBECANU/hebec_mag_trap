%%%%%%% subroutine for calculating magnetic field components at (Ax,Ay,Az)
% for one 'leaf' of a cloverleaf configuration with one loop.

function[Btotal]=extras(k,Rxi,Ryi,Rzi,Rxo,Ryo,Rzo,dthetain,dthetaout,thetain,...
   thetaout,Ax,Ay,Az) 

Binarc=magarcmat(k,Rxi,Ryi,Rzi,Ax,Ay,Az,dthetain,thetain);
Boutarc=magarcmat(k,Rxo,Ryo,Rzo,Ax,Ay,Az,dthetaout,thetaout);

% now lets add all the fields up
Btotal=Binarc+Boutarc;

