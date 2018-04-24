%%%%%%% subroutine for calculating magnetic field components at (Ax,Ay,Az)
% calculates for a line 

function[Bmat]=linemag(k,xstart,xstop,ystart,ystop,zstart,zstop,Axx,Ayy,Azz,intsteps); 

Bmat=zeros(intsteps,3);
xstop=xstop+1e-12;
ystop=ystop+1e-12;
zstop=zstop+1e-12;
[startx,starty,startz,stopx,stopy,stopz]=makelinevector(xstart,xstop,ystart,ystop,zstart,zstop,intsteps) ;

 Bmat=maglinemat(k,startx,starty,startz,stopx,stopy,stopz,Axx,Ayy,Azz); 
  
  Bmat=sum(Bmat,1);











