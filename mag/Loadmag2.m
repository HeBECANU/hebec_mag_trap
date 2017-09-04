%%%%%%% subroutine for calculating magnetic field components at (Ax,Ay,Az)
% for one 'leaf' of a cloverleaf configuration with one loop.

function[Bmat]=Loadmag2(k,x1start,y1start,z1start,x1stop,y1stop,z1stop,...
        x2start,y2start,z2start,x2stop,y2stop,z2stop,x3start,y3start,z3start,x3stop,y3stop,z3stop,...
        x4start,y4start,z4start,x4stop,y4stop,z4stop,Ax,Ay,Az,intsteps)
    
   Bmat=zeros(intsteps,3);
y1stop=y1stop+1e-12;
y2stop=y2stop+1e-12;
z1stop=z1stop+1e-12;
z2stop=z2stop+1e-12;
y3stop=y3stop+1e-12;
y4stop=y4stop+1e-12;
z3stop=z3stop+1e-12;
z4stop=z4stop+1e-12;

   [startx1,starty1,startz1,stopx1,stopy1,stopz1,...
        startx2,starty2,startz2,stopx2,stopy2,stopz2,startx3,starty3,startz3,stopx3,stopy3,stopz3,...
        startx4,starty4,startz4,stopx4,stopy4,stopz4]=makeloadvector(x1start,y1start,z1start,x1stop,y1stop,z1stop,...
        x2start,y2start,z2start,x2stop,y2stop,z2stop,x3start,y3start,z3start,x3stop,y3stop,z3stop,...
        x4start,y4start,z4start,x4stop,y4stop,z4stop,intsteps);

    
    Bmat=Loadmag(k,startx1,starty1,startz1,stopx1,stopy1,stopz1,...
        startx2,starty2,startz2,stopx2,stopy2,stopz2,startx3,starty3,startz3,stopx3,stopy3,stopz3,...
        startx4,starty4,startz4,stopx4,stopy4,stopz4,Ax,Ay,Az,intsteps);
    
     Bmat=sum(Bmat,1);