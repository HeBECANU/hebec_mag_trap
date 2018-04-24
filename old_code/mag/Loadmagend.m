%%%%%%% subroutine for calculating magnetic field components at (Ax,Ay,Az)

function[Bmat]=Loadmagend(k,x1start,y1start,z1start,x1stop,y1stop,z1stop,...
        x2start,y2start,z2start,x2stop,y2stop,z2stop,Ax,Ay,Az,intsteps)
    
   Bmat=zeros(intsteps,3);
y1stop=y1stop+1e-12;
y2stop=y2stop+1e-12;
z1stop=z1stop+1e-12;
z2stop=z2stop+1e-12;

   [startx1,starty1,startz1,stopx1,stopy1,stopz1,...
        startx2,starty2,startz2,stopx2,stopy2,stopz2]=makeloadendvector(x1start,y1start,z1start,x1stop,y1stop,z1stop,...
        x2start,y2start,z2start,x2stop,y2stop,z2stop,intsteps);

    
    Bmat=Loadmagend2(k,startx1,starty1,startz1,stopx1,stopy1,stopz1,...
        startx2,starty2,startz2,stopx2,stopy2,stopz2,Ax,Ay,Az,intsteps);
    
     Bmat=sum(Bmat,1);