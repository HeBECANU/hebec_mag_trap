%THESE ARE THE PRESENT MAGNETIC TRAP COILS at RICE
% the new quads are wrapped to produce zero bias!
% magnetic field calculation program
% calculates the magnetic field using arcs and line segments


% uses two subroutines magarc.m - to calculate the magnetic field components
% for an arc of radius R ...
% magline.m - calculate magnetic field components for  line charge

clear all
%close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%Constants%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
Mu0 = 4*pi*1e-7;
k=(Mu0)/(4*pi);
intsteps=200;

dd=.15;
z1=2*(2.3+dd)/100; % seperation in metres
Nclover=10; % remember the n+1 for the extra coil
Iclover=130; %130
Rin=1.7196/100;% inner radius in m
Rout=4.775/100;% outer radius in m

zbias=2*(3.0+dd)/100;
Rbias=4.6228/100;
Ibias=0;%111.0*.3; %correct value =117.5 ;146.5, 114
Nbias=23;

zcurv=2*(3.0+dd)/100;
Rcurv=2.9972/100;
Icurv=0;%*.3; %130 ; 160
Ncurv=22;
%%%% Note 1T = 10^4 G %%%%%%

sym =1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%First CL calc  .%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% over what range? corrdinates in mm
xmax=1*5;
xmin=-1*5;
ymax=1*5;
ymin=-1*5;
zmax=1*5;
zmin=-1*5;
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
       +(Ibias*halfbiasmag(k,Rbias,zbias,Axx,Ayy,Azz,intsteps))...
+(Icurv*Ncurv*curvmag(k,Rcurv,zcurv,Axx,Ayy,Azz,intsteps))+[0, 0, 0]...
+(Icurv*halfcurvmag(k,Rcurv,zcurv,Axx,Ayy,Azz,intsteps))...
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
       +(Ibias*halfbiasmag(k,Rbias,zbias,Axx,Ayy,Azz,intsteps))...
+(Icurv*Ncurv*curvmag(k,Rcurv,zcurv,Axx,Ayy,Azz,intsteps))+[0, 0, 0]...
+(Icurv*halfcurvmag(k,Rcurv,zcurv,Axx,Ayy,Azz,intsteps))... 
   -(Iclover*extramagwierdsym(k,Rin,Rout,z1,Axx,Ayy,Azz,intsteps,sym)); 
   
end
Bmag=plotline(B,steps,2);
Brad=B;
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
       +(Ibias*halfbiasmag(k,Rbias,zbias,Axx,Ayy,Azz,intsteps))...
+(Icurv*Ncurv*curvmag(k,Rcurv,zcurv,Axx,Ayy,Azz,intsteps))+[0, 0, 0]...
+(Icurv*halfcurvmag(k,Rcurv,zcurv,Axx,Ayy,Azz,intsteps))...
   -(Iclover*extramagwierdsym(k,Rin,Rout,z1,Axx,Ayy,Azz,intsteps,sym)); 
   

end   

Bmag=plotline(B,steps,4);
hold off
pause 
Bmagrad= sqrt(sum(Brad(:,(4:6)).^2,2));
Bmagrad=Bmagrad/1e-4;

kb=1.38e-23;
f = 67e-6; % conversion from Guass to kelvin in |2,2> increases as muB
mass =7/6e23/1000;
Bmagrad=Bmagrad-min(Bmagrad);
U=f*Bmagrad((steps/2)+1:length(Brad(:,2)))*kb;

R=Brad(((steps/2)+1:length(Brad(:,2))),2);
fg(1)=400 ;
fg(2)=1e-7 ;
options(2)=1e-12;
fr=fmins('harmonic',fg,options,[],U,R);

fr
Ufit=(0.5*mass*(2*pi*fr(1))^2*(R-fr(2)).^2);

plot(Brad(:,2)*100,Bmagrad,'b',(R)*100,(Ufit/f/kb),'r')


toc