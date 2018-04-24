%THESE ARE THE PRESENT COILS at RICE
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
tic
Mu0 = 4*pi*1e-7;
k=(Mu0)/(4*pi);
intsteps=200;


zbias=2*2.15/100;
Rbias=3.3/100;%4.6228
Ibias=10; %116.8
Nbias=30;%23.5

zcurv=2*3/100;
Rcurv=2.9972/100;
Icurv=18; %130
Ncurv=22.5;
%%%% Note 1T = 10^4 G %%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%First CL calc  .%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% over what range? corrdinates in mm
xmax=1*10;
xmin=-1*10;
ymax=1*10;
ymin=-1*10;
zmax=1*10;
zmin=-1*10;
steps=20;
[xvec,yvec,zvec]=setupcoords(xmin,xmax,ymin,ymax,zmin,zmax,steps);
counter=0;
B=zeros(21,6);
% first calculate for z-axis
for ttt=1:steps+1
   
   % assign test points
   Axx=0;
   Ayy=0;
   Azz=zvec(ttt);
   
 counter=counter+1;
 B(counter,1)=Axx;
 B(counter,2)=Ayy;
 B(counter,3)=Azz;
 B(counter,(4:6))=(Ibias*Nbias*coilmag(k,Rbias,zbias,Axx,Ayy,Azz,intsteps))...
    -(Ibias*Nbias*coilmag(k,Rbias,-zbias,Axx,Ayy,Azz,intsteps))+[0, 0, 0];   

end

Bmag=plotline(B,steps,3);

hold on
counter=0;
% now calculate for x-axis
for ttt=1:steps+1
   
   % assign test points
   Axx=0;
   Ayy=yvec(ttt);
   Azz=0;
   
 counter=counter+1;
 B(counter,1)=Axx;
 B(counter,2)=Ayy;
 B(counter,3)=Azz;
 B(counter,(4:6))=(Ibias*Nbias*coilmag(k,Rbias,zbias,Axx,Ayy,Azz,intsteps))...
    -(Ibias*Nbias*coilmag(k,Rbias,-zbias,Axx,Ayy,Azz,intsteps))+[0, 0, 0];   
end
Bmag=plotline(B,steps,2);
counter=0;

% now calculate for xy-axis
for ttt=1:steps+1   
   % assign test points
   Axx=xvec(ttt);
   Ayy=yvec(ttt);
   Azz=0;
   
 counter=counter+1;
 B(counter,1)=Axx;
 B(counter,2)=Ayy;
 B(counter,3)=Azz;
 B(counter,(4:6))=(Ibias*Nbias*coilmag(k,Rbias,zbias,Axx,Ayy,Azz,intsteps))...
    -(Ibias*Nbias*coilmag(k,Rbias,-zbias,Axx,Ayy,Azz,intsteps))+[0, 0, 0];   
end   

Bmag=plotline(B,steps,4);
hold off
toc