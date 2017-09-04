%THESE ARE THE PRESENT MOT COILS at ANU
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
k=(Mu0)/(4*pi);
intsteps=400;

%%%%%%%%%%%%%%%%%%%%%%%For coils OD =1.8mm, ID =0.7 mm%%%%%%%%%%%%%
%%%%% R= 21/1000 Ohms/m
I1=10; %

xstart=-.5/100; xstop=.5/100;
zstart=-1/100; zstop=0/100;
ystart=-2/1000; ystop=-2/1000;
bias=[ -2/1e4 0 0/1e4]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%First CL calc  .%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% over what range? corrdinates in mm
range=4;
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

tp=-.0/100;

% first calculate for z-axis
for ttt=1:steps+1
   
   % assign test points
   Axx=.0;
   Ayy=.0;
   Azz=zvec(ttt);
   
 counter=counter+1;
 B(counter,1)=Axx;
 B(counter,2)=Ayy;
 B(counter,3)=Azz;

B(counter,(4:6))=B(counter,(4:6))+(I1*linemag(k,xstart,xstart,ystart,ystart,zstart,zstop,Axx,Ayy,Azz,intsteps)); 
B(counter,(4:6))=B(counter,(4:6))+(I1*linemag(k,xstart,xstop,ystart,ystart,zstop,zstop,Axx,Ayy,Azz,intsteps)); 
B(counter,(4:6))=B(counter,(4:6))+(I1*linemag(k,xstop,xstop,ystart,ystart,zstop,zstart,Axx,Ayy,Azz,intsteps))+bias; 

end

Bmag=plotline(B,steps,3);
hold on
pause

counter=0;
% now calculate for y-axis
for ttt=1:steps+1
   
   % assign test points
   Axx=.0;
   Ayy=yvec(ttt);
   Azz=0;
   
 counter=counter+1;
 B(counter,1)=Axx;
 B(counter,2)=Ayy;
 B(counter,3)=Azz;
B(counter,(4:6))=B(counter,(4:6))+(I1*linemag(k,xstart,xstart,ystart,ystart,zstart,zstop,Axx,Ayy,Azz,intsteps)); 
B(counter,(4:6))=B(counter,(4:6))+(I1*linemag(k,xstart,xstop,ystart,ystart,zstop,zstop,Axx,Ayy,Azz,intsteps)); 
B(counter,(4:6))=B(counter,(4:6))+(I1*linemag(k,xstop,xstop,ystart,ystart,zstop,zstart,Axx,Ayy,Azz,intsteps))+bias; 

end
Bmag=plotline(B,steps,2);
counter=0;


% now calculate for x-axis
for ttt=1:steps+1
   
   % assign test points
   Axx=xvec(ttt);
   Ayy=.0;
   Azz=0;
   
 counter=counter+1;
 B(counter,1)=Axx;
 B(counter,2)=Ayy;
 B(counter,3)=Azz;
B(counter,(4:6))=B(counter,(4:6))+(I1*linemag(k,xstart,xstart,ystart,ystart,zstart,zstop,Axx,Ayy,Azz,intsteps)); 
B(counter,(4:6))=B(counter,(4:6))+(I1*linemag(k,xstart,xstop,ystart,ystart,zstop,zstop,Axx,Ayy,Azz,intsteps)); 
B(counter,(4:6))=B(counter,(4:6))+(I1*linemag(k,xstop,xstop,ystart,ystart,zstop,zstart,Axx,Ayy,Azz,intsteps))+bias; 
 



end

%Bmag=plotline(B,steps,1);
counter=0;
AA=B;

% now calculate for xy-axis
for ttt=1:steps+1   
   % assign test points
  Axx=xvec(ttt);
   Ayy=yvec(ttt);
   Azz=tp;
   
 counter=counter+1;
 B(counter,1)=Axx;
 B(counter,2)=Ayy;
 B(counter,3)=Azz;
B(counter,(4:6))=B(counter,(4:6))+(I1*linemag(k,xstart,xstart,ystart,ystart,zstart,zstop,Axx,Ayy,Azz,intsteps)); 
B(counter,(4:6))=B(counter,(4:6))+(I1*linemag(k,xstart,xstop,ystart,ystart,zstop,zstop,Axx,Ayy,Azz,intsteps)); 
B(counter,(4:6))=B(counter,(4:6))+(I1*linemag(k,xstop,xstop,ystart,ystart,zstop,zstart,Axx,Ayy,Azz,intsteps))+bias; 


end   

Bmag=plotline(B,steps,4);
key
grid on
%hold off
%plot(x(:,1),x(:,5),'kh',xy(:,1),xy(:,5),'yo',z(:,1),z(:,5),'c.',y(:,1),y(:,5),'rs')
hold off

toc