%%%%%%% subroutine for calculating magnetic field components at (Ax,Ay,Az)
% for one 'leaf' of a cloverleaf configuration with one loop.

function[xvec,yvec,zvec]=setupcoords(xmin,xmax,ymin,ymax,zmin,zmax,steps); 

%the actual step size in metres
xstep=(xmax-xmin)/(1000*steps);
ystep=(ymax-ymin)/(1000*steps);
zstep=(zmax-zmin)/(1000*steps);

%coordinates of test points
xvec =xmin/1000:xstep:xmax/1000;
yvec =ymin/1000:ystep:ymax/1000;
zvec =zmin/1000:zstep:zmax/1000;




