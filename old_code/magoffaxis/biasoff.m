%%%%%%% subroutine for calculating magnetic field components at (Ax,Ay,Az)
% for one 'leaf' of a cloverleaf configuration with one loop.

function[Btotal]=biasoff(k,offset,Rx,Ry,Rz,dtheta,theta,Axx,Ayy,Azz) 

Btotal=magarcmatoffaxis(k,offset,Rx,Ry,Rz,Axx,Ayy,Azz,dtheta,theta);

