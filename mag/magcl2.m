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

z1=.28; % seperation in metres
intsteps=4;
Rin=3/1000;% inner radius in mm
Rout=4/1000;% outer radius in mm
dtheta=(pi/2)/intsteps;

%%%% Note 1T = 10^4 G %%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%First CL calc  .%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%first define clover leaf coordinates
% define line coordinates first
x1start=-Rin;
x1stop=-Rout;
y1start=0;
y1stop=1e-20;
z1start=-z1/2;
z1stop=-z1/2+1e-12;
dx1=(x1stop-x1start)/intsteps;
dy1=(y1stop-y1start)/intsteps;
dz1=(z1stop-z1start)/intsteps;
startx1=x1start:dx1:x1stop-dx1;
starty1=y1start:dy1:y1stop-dy1;
startz1=z1start:dz1:z1stop-dz1;
stopx1=x1start+dx1:dx1:x1stop;
stopy1=y1start+dy1:dy1:y1stop;
stopz1=z1start+dz1:dz1:z1stop;

x2start=0;
x2stop=1e-20;
y2start=-Rout;
y2stop=-Rin;
z2start=-z1/2;
z2stop=-z1/2+1e-12;
dx2=(x2stop-x2start)/intsteps;
dy2=(y2stop-y2start)/intsteps;
dz2=(z2stop-z2start)/intsteps;
startx2=x2start:dx2:x2stop-dx2;
starty2=y2start:dy2:y2stop-dy2;
startz2=z2start:dz2:z2stop-dz2;
stopx2=x2start+dx2:dx2:x2stop;
stopy2=y2start+dy2:dy2:y2stop;
stopz2=z2start+dz2:dz2:z2stop;

% now arc coordinates
% remember to get the coordinates right in terms of direction of the current
theta1 = (pi)+dtheta/2:dtheta:(3*pi/2)-dtheta/2;
   
   Rxi=Rin*cos(theta1);
   Ryi=Rin*sin(theta1);
   Rzi=-z1/2;
   
   theta2 = (3*pi/2)-dtheta/2:-dtheta:(pi)+dtheta/2;
   
   Rxo=Rout*cos(theta2);
   Ryo=Rout*sin(theta2);
   Rzo=-z1/2;

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

counter=0;
% Bloop=k*2*pi*R1^2*(Rz^2+R1^2)^(-1.5)
for ttt=1:zzstep
   Bmat=zeros(intsteps,3);

   
   Bmat=clmagmat(k,startx1,starty1,startz1,stopx1,stopy1,...
   stopz1,startx2,starty2,startz2,stopx2,stopy2,stopz2,Rxi,Ryi,...
   Rzi,Rxo,Ryo,Rzo,dtheta,theta1,theta2,Ax,Ay,Az(ttt));

 Bmat=sum(Bmat,1);
 counter=counter+1;
 B(counter,1)=Ax;
 B(counter,2)=Ay;
 B(counter,3)=Az(ttt);
 B(counter,(4:6))=Bmat;
end

Bmag=plotline(B,zzstep);
