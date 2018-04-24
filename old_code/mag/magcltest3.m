% magnetic field calculation program
% calculates the magnetic field using arcs and line segments


% uses two subroutines magarc.m - to calculate the magnetic field components
% for an arc of radius R ...
% magline.m - calculate magnetic field components for  line charge

clear all
% close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%Constants%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
Mu0 = 4*pi*1e-7;
I=1; % in Amps note positive current indicates anti clockise current flow
	  % when viewed from the origin
N=1;
k=(Mu0*N*I)/(4*pi);

z1=2.3/100; % seperation in metres
intsteps=200;
Rin=1.7196/100;% inner radius in m
Rout=4.775/100;% outer radius in m
dtheta=(pi/2)/intsteps;

%%%% Note 1T = 10^4 G %%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%First CL calc  .%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%First define region of interest%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%How many steps?

xstep=100;
zstep=100;
ystep=100;
zzstep=zstep;
% over what range? corrdinates in mm
xmax=20;
xmin=-20;
ymax=10;
ymin=-10;
zmax=50;
zmin=-50;

%the actual step size in metres
xstep=(xmax-xmin)/(1000*xstep);
ystep=(ymax-ymin)/(1000*ystep);
zstep=(zmax-zmin)/(1000*zstep);

%coordinates of test points
Ax =xmin/1000:xstep:xmax/1000;
Ay =0;%ymin/1000:ystep:ymax/1000;
Az =0;%zmin/1000:zstep:zmax/1000;

counter=0;

for ttt=1:zzstep
   Bmat=zeros(intsteps,3);
   
   % assign test points
   Axx=Ax(ttt);
   Ayy=Ay;
   Azz=Az;
   
   % cloverleaf 1
  x1start=-Rout;,x1stop=-Rin;,y1start=0;,y1stop=1e-20;,z1start=-z1/2;,z1stop=-z1/2+1e-12;
x2start=0;,x2stop=1e-20;,y2start=Rin;,y2stop=Rout;,z2start=-z1/2;,z2stop=-z1/2+1e-12;
thetainstart=pi;,thetainstop=pi/2;,thetaoutstart=pi/2;,thetaoutstop=pi;
[startx1,starty1,startz1,stopx1,stopy1,stopz1,startx2,starty2,startz2,stopx2,stopy2,...
      stopz2,Rxi,Ryi,Rzi,Rxo,Ryo,Rzo,dtheta1,dtheta2,theta1,theta2]=makevector(x1start,...
   x1stop,y1start,y1stop,z1start,z1stop,x2start,x2stop,y2start,y2stop,z2start,z2stop,...
Rin,Rout,thetaoutstart,thetaoutstop,thetainstart,thetainstop,intsteps); 

  Bmat=clmagmat(k,startx1,starty1,startz1,stopx1,stopy1,...
   stopz1,startx2,starty2,startz2,stopx2,stopy2,stopz2,Rxi,Ryi,...
   Rzi,Rxo,Ryo,Rzo,dtheta1,dtheta2,theta1,theta2,Axx,Ayy,Azz);

x1start=-Rout;,x1stop=-Rin;,y1start=0;,y1stop=1e-20;,z1start=-z1/2;,z1stop=-z1/2+1e-12;
x2start=0;,x2stop=1e-20;,y2start=-Rin;,y2stop=-Rout;,z2start=-z1/2;,z2stop=-z1/2+1e-12;
thetainstart=pi;,thetainstop=3*pi/2;,thetaoutstart=3*pi/2;,thetaoutstop=pi;
[startx1,starty1,startz1,stopx1,stopy1,stopz1,startx2,starty2,startz2,stopx2,stopy2,...
      stopz2,Rxi,Ryi,Rzi,Rxo,Ryo,Rzo,dtheta1,dtheta2,theta1,theta2]=makevector(x1start,...
   x1stop,y1start,y1stop,z1start,z1stop,x2start,x2stop,y2start,y2stop,z2start,z2stop,...
   Rin,Rout,thetaoutstart,thetaoutstop,thetainstart,thetainstop,intsteps); 

  Bmat=Bmat+clmagmat(k,startx1,starty1,startz1,stopx1,stopy1,...
   stopz1,startx2,starty2,startz2,stopx2,stopy2,stopz2,Rxi,Ryi,...
   Rzi,Rxo,Ryo,Rzo,dtheta1,dtheta2,theta1,theta2,Axx,Ayy,Azz);

 x1start=Rout;,x1stop=Rin;,y1start=0;,y1stop=1e-20;,z1start=-z1/2;,z1stop=-z1/2+1e-12;
x2start=0;,x2stop=1e-20;,y2start=Rin;,y2stop=Rout;,z2start=-z1/2;,z2stop=-z1/2+1e-12;
thetainstart=0;,thetainstop=pi/2;,thetaoutstart=pi/2;,thetaoutstop=0;
[startx1,starty1,startz1,stopx1,stopy1,stopz1,startx2,starty2,startz2,stopx2,stopy2,...
      stopz2,Rxi,Ryi,Rzi,Rxo,Ryo,Rzo,dtheta1,dtheta2,theta1,theta2]=makevector(x1start,...
   x1stop,y1start,y1stop,z1start,z1stop,x2start,x2stop,y2start,y2stop,z2start,z2stop,...
Rin,Rout,thetaoutstart,thetaoutstop,thetainstart,thetainstop,intsteps); 

   Bmat=Bmat+clmagmat(k,startx1,starty1,startz1,stopx1,stopy1,...
   stopz1,startx2,starty2,startz2,stopx2,stopy2,stopz2,Rxi,Ryi,...
   Rzi,Rxo,Ryo,Rzo,dtheta1,dtheta2,theta1,theta2,Axx,Ayy,Azz);


x1start=Rout;,x1stop=Rin;,y1start=0;,y1stop=1e-20;,z1start=-z1/2;,z1stop=-z1/2+1e-12;
x2start=0;,x2stop=1e-20;,y2start=-Rin;,y2stop=-Rout;,z2start=-z1/2;,z2stop=-z1/2+1e-12;
thetainstart=2*pi;,thetainstop=3*pi/2;,thetaoutstart=3*pi/2;,thetaoutstop=2*pi;
[startx1,starty1,startz1,stopx1,stopy1,stopz1,startx2,starty2,startz2,stopx2,stopy2,...
      stopz2,Rxi,Ryi,Rzi,Rxo,Ryo,Rzo,dtheta1,dtheta2,theta1,theta2]=makevector(x1start,...
   x1stop,y1start,y1stop,z1start,z1stop,x2start,x2stop,y2start,y2stop,z2start,z2stop,...
Rin,Rout,thetaoutstart,thetaoutstop,thetainstart,thetainstop,intsteps); 

  Bmat=Bmat+clmagmat(k,startx1,starty1,startz1,stopx1,stopy1,...
 stopz1,startx2,starty2,startz2,stopx2,stopy2,stopz2,Rxi,Ryi,...
 Rzi,Rxo,Ryo,Rzo,dtheta1,dtheta2,theta1,theta2,Axx,Ayy,Azz);

 % cloverleaf 2
  x1start=-Rin;,x1stop=-Rout;,y1start=0;,y1stop=1e-20;,z1start=z1/2;,z1stop=z1/2+1e-12;
x2start=0;,x2stop=1e-20;,y2start=Rout;,y2stop=Rin;,z2start=z1/2;,z2stop=z1/2+1e-12;
thetainstart=pi/2;,thetainstop=pi;,thetaoutstart=pi;,thetaoutstop=pi/2;
[startx1,starty1,startz1,stopx1,stopy1,stopz1,startx2,starty2,startz2,stopx2,stopy2,...
      stopz2,Rxi,Ryi,Rzi,Rxo,Ryo,Rzo,dtheta1,dtheta2,theta1,theta2]=makevector(x1start,...
   x1stop,y1start,y1stop,z1start,z1stop,x2start,x2stop,y2start,y2stop,z2start,z2stop,...
Rin,Rout,thetaoutstart,thetaoutstop,thetainstart,thetainstop,intsteps); 

  Bmat=Bmat+clmagmat(k,startx1,starty1,startz1,stopx1,stopy1,...
   stopz1,startx2,starty2,startz2,stopx2,stopy2,stopz2,Rxi,Ryi,...
   Rzi,Rxo,Ryo,Rzo,dtheta1,dtheta2,theta1,theta2,Axx,Ayy,Azz);

x1start=-Rin;,x1stop=-Rout;,y1start=0;,y1stop=1e-20;,z1start=z1/2;,z1stop=z1/2+1e-12;
x2start=0;,x2stop=1e-20;,y2start=-Rout;,y2stop=-Rin;,z2start=z1/2;,z2stop=z1/2+1e-12;
thetainstart=3*pi/2;,thetainstop=pi;,thetaoutstart=pi;,thetaoutstop=3*pi/2;
[startx1,starty1,startz1,stopx1,stopy1,stopz1,startx2,starty2,startz2,stopx2,stopy2,...
      stopz2,Rxi,Ryi,Rzi,Rxo,Ryo,Rzo,dtheta1,dtheta2,theta1,theta2]=makevector(x1start,...
   x1stop,y1start,y1stop,z1start,z1stop,x2start,x2stop,y2start,y2stop,z2start,z2stop,...
   Rin,Rout,thetaoutstart,thetaoutstop,thetainstart,thetainstop,intsteps); 

  Bmat=Bmat+clmagmat(k,startx1,starty1,startz1,stopx1,stopy1,...
   stopz1,startx2,starty2,startz2,stopx2,stopy2,stopz2,Rxi,Ryi,...
   Rzi,Rxo,Ryo,Rzo,dtheta1,dtheta2,theta1,theta2,Axx,Ayy,Azz);

 x1start=Rin;,x1stop=Rout;,y1start=0;,y1stop=1e-20;,z1start=z1/2;,z1stop=z1/2+1e-12;
x2start=0;,x2stop=1e-20;,y2start=Rout;,y2stop=Rin;,z2start=z1/2;,z2stop=z1/2+1e-12;
thetainstart=pi/2;,thetainstop=0;,thetaoutstart=0;,thetaoutstop=pi/2;
[startx1,starty1,startz1,stopx1,stopy1,stopz1,startx2,starty2,startz2,stopx2,stopy2,...
      stopz2,Rxi,Ryi,Rzi,Rxo,Ryo,Rzo,dtheta1,dtheta2,theta1,theta2]=makevector(x1start,...
   x1stop,y1start,y1stop,z1start,z1stop,x2start,x2stop,y2start,y2stop,z2start,z2stop,...
Rin,Rout,thetaoutstart,thetaoutstop,thetainstart,thetainstop,intsteps); 

   Bmat=Bmat+clmagmat(k,startx1,starty1,startz1,stopx1,stopy1,...
   stopz1,startx2,starty2,startz2,stopx2,stopy2,stopz2,Rxi,Ryi,...
   Rzi,Rxo,Ryo,Rzo,dtheta1,dtheta2,theta1,theta2,Axx,Ayy,Azz);


x1start=Rin;,x1stop=Rout;,y1start=0;,y1stop=1e-20;,z1start=z1/2;,z1stop=z1/2+1e-12;
x2start=0;,x2stop=1e-20;,y2start=-Rout;,y2stop=-Rin;,z2start=z1/2;,z2stop=z1/2+1e-12;
thetainstart=3*pi/2;,thetainstop=2*pi;,thetaoutstart=2*pi;,thetaoutstop=3*pi/2;
[startx1,starty1,startz1,stopx1,stopy1,stopz1,startx2,starty2,startz2,stopx2,stopy2,...
      stopz2,Rxi,Ryi,Rzi,Rxo,Ryo,Rzo,dtheta1,dtheta2,theta1,theta2]=makevector(x1start,...
   x1stop,y1start,y1stop,z1start,z1stop,x2start,x2stop,y2start,y2stop,z2start,z2stop,...
Rin,Rout,thetaoutstart,thetaoutstop,thetainstart,thetainstop,intsteps); 

  Bmat=Bmat+clmagmat(k,startx1,starty1,startz1,stopx1,stopy1,...
 stopz1,startx2,starty2,startz2,stopx2,stopy2,stopz2,Rxi,Ryi,...
 Rzi,Rxo,Ryo,Rzo,dtheta1,dtheta2,theta1,theta2,Axx,Ayy,Azz);




Bmat=sum(Bmat,1);
 counter=counter+1;
 B(counter,1)=Axx;
 B(counter,2)=Ayy;
 B(counter,3)=Azz;
 B(counter,(4:6))=Bmat;
end
Bmat;
Bmag=plotline(B,zzstep);
toc