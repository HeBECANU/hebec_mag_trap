%%%%%%% subroutine for making the relevant cloverleaf coordinates
function[startx,starty,startz,stopx,stopy,stopz]=makelinevector(xstart,xstop,ystart,ystop,zstart,zstop,intsteps) 

% define line coordinates 
dx=(xstop-xstart)/intsteps;
dy=(ystop-ystart)/intsteps;
dz=(zstop-zstart)/intsteps;
startx=xstart:dx:xstop-dx;
starty=ystart:dy:ystop-dy;
startz=zstart:dz:zstop-dz;
stopx=xstart+dx:dx:xstop;
stopy=ystart+dy:dy:ystop;
stopz=zstart+dz:dz:zstop;




