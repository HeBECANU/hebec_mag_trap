%THESE ARE THE PRESENT MAGNETIC TRAP COILS at RICE
% before the quad's blew up and were replaced
% the new quads are wrapped to produce zero bias!
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


z1=2*2.3/100; % seperation in metres
Nclover=14; % remember the n+1 for the extra coil
Iclover=200; %130
Rin=1.7196/100;% inner radius in m
Rout=4.775/100;% outer radius in m
sym=1; %0 % asymmetry in extra clover wraps

zbias=2*3.0/100;
Rbias=6/100;%4.6228/100;
Ibias=0;%118.5; %correct value =118.5 ;146.5 114.5
Nbias=23;%23.5

zcurv=2*3.0/100;
Rcurv=2.9972/100;
Icurv=130; %130 ; 160
Ncurv=22;%22.5
%%%% Note 1T = 10^4 G %%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%First CL calc  .%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% over what range? corrdinates in mm
range=100;
xmax=1*range;
xmin=-1*range;
ymax=1*range;
ymin=-1*range;
zmax=1*range;
zmin=-1*range;
steps=200;
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
 B(counter,(4:6))= (Iclover*Nclover*clovermag(k,Rin,Rout,z1,Axx,Ayy,Azz,intsteps,0,0))...
    +(Ibias*Nbias*biasmag(k,Rbias,zbias,Axx,Ayy,Azz,intsteps))...
    +(Icurv*Ncurv*curvmag(k,Rcurv,zcurv,Axx,Ayy,Azz,intsteps))+[0, 0, 0]...
   -(Iclover*extramagwierdsym(k,Rin,Rout,z1,Axx,Ayy,Azz,intsteps,sym)); 
   
   

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
 B(counter,(4:6))=(Iclover*Nclover*clovermag(k,Rin,Rout,z1,Axx,Ayy,Azz,intsteps,0,0))...
	+(Ibias*Nbias*biasmag(k,Rbias,zbias,Axx,Ayy,Azz,intsteps))...
    +(Icurv*Ncurv*curvmag(k,Rcurv,zcurv,Axx,Ayy,Azz,intsteps))+[0, 0, 0]... 
   -(Iclover*extramagwierdsym(k,Rin,Rout,z1,Axx,Ayy,Azz,intsteps,sym)); 
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
 B(counter,(4:6))=(Iclover*Nclover*clovermag(k,Rin,Rout,z1,Axx,Ayy,Azz,intsteps,0,0))...
	+(Ibias*Nbias*biasmag(k,Rbias,zbias,Axx,Ayy,Azz,intsteps))...
    +(Icurv*Ncurv*curvmag(k,Rcurv,zcurv,Axx,Ayy,Azz,intsteps))+[0, 0, 0]... 
   -(Iclover*extramagwierdsym(k,Rin,Rout,z1,Axx,Ayy,Azz,intsteps,sym)); 
end   

Bmag=plotline(B,steps,4);
hold off
toc