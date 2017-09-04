% magnetic field calculation program
% calculates the magnetic field using arcs and line segments


% uses two subroutines magarc.m - to calculate the magnetic field components
% for an arc of radius R ...
% magline.m - calculate magnetic field components for  line charge

clear all
close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%Constants%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Mu0 = 4*pi*1e-7;
I=1; % in Amps note positive current indicates anti clockise current flow
	  % when viewed from the origin
N=1;
k=(Mu0*N*I)/(4*pi);
R1 =.14; % in metres
z1=.28; % seperation in metres
intsteps=200;

%%%% Note 1T = 10^4 G %%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%First line calc.%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%first define line coordinates

xlinestart=-.1;
xlinestop=.1% in metres
dx=(xlinestop-xlinestart)/intsteps;
ylinestart=0;
ylinestop=1e-20;
dy=(ylinestop-ylinestart)/intsteps;
zlinestart=0;
zlinestop=1e-20;
dz=(zlinestop-zlinestart)/intsteps;
startx=xlinestart:dx:xlinestop-dx;
starty=ylinestart:dy:ylinestop-dy;
startz=zlinestart:dz:zlinestop-dz;
stopx=xlinestart+dx:dx:xlinestop;
stopy=ylinestart+dy:dy:ylinestop;
stopz=zlinestart+dz:dz:zlinestop;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%First define region of interest%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%How many steps?

ystep=50;
zstep=100;
yystep=ystep;
zzstep=zstep;
numpoints=(yystep)*(zzstep);
% over what range? corrdinates in mm
xmax=1;
xmin=-1;
ymax=10;
ymin=-10;
zmax=50;
zmin=-50;

%the actual step size in metres
% xstep=(xmax-xmin)/(1000*xstep);
%ystep=(ymax-ymin)/(1000*ystep);
zstep=(zmax-zmin)/(1000*zstep);

%coordinates of test points
Ax =0;%xmin/1000:xstep:xmax/1000;
Ay =0;%ymin/1000:ystep:ymax/1000;
Az =zmin/1000:zstep:zmax/1000;

Bmat=maglinemat(k,startx,starty,startz,stopx,stopy,stopz,0,0,.1);
BB=sum(Bmat,1);
BBB=sqrt(sum(BB.^2,2))
Bline=k*sqrt(2)/.1


tic
counter=0;
% Bloop=k*2*pi*R1^2*(Rz^2+R1^2)^(-1.5)
for ttt=1:zzstep
   Bmat=zeros(intsteps,3);

   
   Bmat=maglinemat(k,startx,starty,startz,stopx,stopy,stopz,Ax,Ay,Az(ttt));
 
 Bmat=sum(Bmat,1);
 counter=counter+1;
 B(counter,1)=Ax;
 B(counter,2)=Ay;
 B(counter,3)=Az(ttt);
 B(counter,(4:6))=Bmat;
end
toc

Bmag=plotline(B,zzstep);
