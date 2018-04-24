%%%%%%% subroutine for calculating magnetic field components at (Ax,Ay,Az)
% for one 'leaf' of a cloverleaf configuration with one loop.

function[Btotal]=clmagmatextra(k,startx1,starty1,startz1,stopx1,stopy1,...
   stopz1,startx2,starty2,startz2,stopx2,stopy2,stopz2,Rxi,Ryi,...
   Rzi,Rxo,Ryo,Rzo,dthetain,dthetaout,thetain,thetaout,Ax,Ay,Az) 

Boutarc=magarcmat(k,Rxo,Ryo,Rzo,Ax,Ay,Az,dthetaout,thetaout);

Btotal=Boutarc;

