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
R1 =.025; % in metres
z1=.025; % seperation in metres
intsteps=200;
dtheta=(2*pi)/intsteps;

%%%% Note 1T = 10^4 G %%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%First coil calc.%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%First define region of interest%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%How many steps?

ystep=20;
zstep=20;
yystep=ystep+1;
zzstep=zstep+1;
numpoints=(yystep)*(zzstep);
% over what range? corrdinates in mm
xmax=1;
xmin=-1;
ymax=50;
ymin=-50;
zmax=50;
zmin=-50;

%the actual step size in metres
% xstep=(xmax-xmin)/(1000*xstep);
ystep=(ymax-ymin)/(1000*ystep);
zstep=(zmax-zmin)/(1000*zstep);

%coordinates of test points
Ax =0;%xmin/1000:xstep:xmax/1000;
Ay =ymin/1000:ystep:ymax/1000;
Az =zmin/1000:zstep:zmax/1000;

tic
counter=0;
% Bloop=k*2*pi*R1^2*(Rz^2+R1^2)^(-1.5)
for ttt=1:zzstep
for tt=1:yystep
   Bmat=zeros(intsteps,3);

theta = dtheta/2:dtheta:(2*pi)-dtheta/2;
   
   Rx=R1*cos(theta);
   Ry=R1*sin(theta);
   Rz=-z1/2;
   
   Bmat=magarcmat(k,Rx,Ry,Rz,Ax,Ay(tt),Az(ttt),dtheta,theta);
 
 Bmat=sum(Bmat,1);
 counter=counter+1;
 B(counter,1)=Ax;
  B(counter,2)=Ay(tt);
 B(counter,3)=Az(ttt);
 B(counter,(4:6))=Bmat;
end
end
 %B1=B(:,(4:6));

toc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%Second coil calculation%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%If you don't want to use any symmetry put dumb =1
%%%%% i.e. calculate for each coil
dumb=0;
if dumb ==1
counter=0;
 for ttt=1:zzstep
for tt=1:yystep
   Bmat=zeros(intsteps,3);

theta = dtheta/2:dtheta:(2*pi)-dtheta/2;
   
   Rx=R1*cos(theta);
   Ry=R1*sin(theta);
   Rz=z1/2;
   
   Bmat=magarcmat(k,Rx,Ry,Rz,Ax,Ay(tt),Az(ttt),dtheta,theta);
 
 Bmat=sum(Bmat,1);
 counter=counter+1;
 B(counter,(4:6))=B(counter,(4:6))+Bmat;
end
end
elseif dumb==0
   B2=flipud(B(:,(4:6)));
   B(:,(4:6))=B(:,(4:6))+B2;
   end
 	plotmag3(B,yystep,zzstep);