%%%%%%% subroutine for calculating magnetic field components at (Ax,Ay,Az)

function[Btotal]=extraswierd(k,Rxi,Ryi,Rzi,Rxo,Ryo,Rzo,dthetain,dthetaout,thetain,...
   thetaout,Ax,Ay,Az) 

Binarc=0;
Boutarc=magarcmat(k,Rxo,Ryo,Rzo,Ax,Ay,Az,dthetaout,thetaout);

Btotal=Boutarc;

