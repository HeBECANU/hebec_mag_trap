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
     Nclover=10;
     Iclover=260;
k=(Mu0)/(4*pi);

z1=2*2.3/100; % seperation in metres
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
xmax=10;
xmin=-10;
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
   
   % assign test points
   Axx=Ax(ttt);
   Ayy=Ay;
   Azz=Az;
   
 counter=counter+1;
 B(counter,1)=Axx;
 B(counter,2)=Ayy;
 B(counter,3)=Azz;
 B(counter,(4:6))=(Iclover*(Nclover+1)*clovermag(k,Rin,Rout,z1,Axx,Ayy,Azz,intsteps,0))...
    +(Iclover*extramagwierd(k,Rin,Rout,z1,Axx,Ayy,Azz,intsteps)); 



end

Bmag=plotline(B,zzstep);
toc