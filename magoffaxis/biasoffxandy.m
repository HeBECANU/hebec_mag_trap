%%%%%%% subroutine for calculating magnetic field components at (Ax,Ay,Az)
% for one 'leaf' of a cloverleaf configuration with one loop.

function[Btotal]=biasoffxandy(k,offsetx,offsety,Rx,Ry,Rz,dtheta,theta,Axx,Ayy,Azz) 

Btotal=magarcmatoffaxisxandy(k,offsetx,offsety,Rx,Ry,Rz,Axx,Ayy,Azz,dtheta,theta);

