%%%%%%% subroutine for calculating magnetic field components at (Ax,Ay,Az)
% for one 'leaf' of a cloverleaf configuration with one loop.

function[Btotal]=biasz(k,Rx,Ry,Rz,dtheta,theta,Axx,Ayy,Azz) 

Btotal=magarcmatx(k,Rx,Ry,Rz,Axx,Ayy,Azz,dtheta,theta);


