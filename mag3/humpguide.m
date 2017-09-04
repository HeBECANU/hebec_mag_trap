%THESE ARE THE PRESENT MOT COILS at ANU
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
intsteps=400;

%%%%%%%%%%%%%%%%%%%%%%%For coils OD =1.8mm, ID =0.7 mm%%%%%%%%%%%%%
%%%%% R= 21/1000 Ohms/m
I=1; %-60


xc=.5/100;
yc=.5/100;
zc=.5/100;
ystart1=yc;xstart1=xc;zstart1=-zc;ystop1=yc;xstop1=xc;zstop1=zc;
ystart2=-yc;xstart2=xc;zstart2=-zc;ystop2=-yc;xstop2=xc;zstop2=zc;
ystart3=-yc;xstart3=-xc;zstart3=-zc;ystop3=-yc;xstop3=-xc;zstop3=zc;
ystart4=yc;xstart4=-xc;zstart4=-zc;ystop4=yc;xstop4=-xc;zstop4=zc;

%%%% Note 1T = 10^4 G %%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%First CL calc  .%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% over what range? corrdinates in mm
range=5;
xmax=1*range;
xmin=-1*range;
ymax=1*range;
ymin=-1*range;
zmax=1*range;
zmin=-1*range;
steps=100;
[xvec,yvec,zvec]=setupcoords(xmin,xmax,ymin,ymax,zmin,zmax,steps);
counter=0;
B=zeros(21,6);

tp=-.0/100;

% first calculate for z-axis
for ttt=1:steps+1
   
   % assign test points
   Axx=tp;
   Ayy=tp;
   Azz=zvec(ttt);
   
 counter=counter+1;
 B(counter,1)=Axx;
 B(counter,2)=Ayy;
 B(counter,3)=Azz;
 B(counter,(4:6))= I*(linemag(k,xstart1,xstop1,ystart1,ystop1,zstart1,zstop1,Axx,Ayy,Azz,intsteps)-...
                        linemag(k,xstart2,xstop2,ystart2,ystop2,zstart2,zstop2,Axx,Ayy,Azz,intsteps)+...
                        linemag(k,xstart3,xstop3,ystart3,ystop3,zstart3,zstop3,Axx,Ayy,Azz,intsteps)-...
                        linemag(k,xstart4,xstop4,ystart4,ystop4,zstart4,zstop4,Axx,Ayy,Azz,intsteps)); 
end
%Bmag=plotline(B,steps,3);
hold on


counter=0;
% now calculate for y-axis
for ttt=1:steps+1
   
   % assign test points
   Axx=tp;
   Ayy=yvec(ttt);
   Azz=0;
   
 counter=counter+1;
 B(counter,1)=Axx;
 B(counter,2)=Ayy;
 B(counter,3)=Azz;
 B(counter,(4:6))= I*(linemag(k,xstart1,xstop1,ystart1,ystop1,zstart1,zstop1,Axx,Ayy,Azz,intsteps)-...
                        linemag(k,xstart2,xstop2,ystart2,ystop2,zstart2,zstop2,Axx,Ayy,Azz,intsteps)+...
                        linemag(k,xstart3,xstop3,ystart3,ystop3,zstart3,zstop3,Axx,Ayy,Azz,intsteps)-...
                        linemag(k,xstart4,xstop4,ystart4,ystop4,zstart4,zstop4,Axx,Ayy,Azz,intsteps)); 

end

Bmag=plotline(B,steps,2);
counter=0;
break
% now calculate for x-axis
for ttt=1:steps+1
   
   % assign test points
   Axx=xvec(ttt);
   Ayy=tp;
   Azz=0;
   
 counter=counter+1;
 B(counter,1)=Axx;
 B(counter,2)=Ayy;
 B(counter,3)=Azz;
 B(counter,(4:6))= I*(linemag(k,xstart1,xstop1,ystart1,ystop1,zstart1,zstop1,Axx,Ayy,Azz,intsteps)-...
                        linemag(k,xstart2,xstop2,ystart2,ystop2,zstart2,zstop2,Axx,Ayy,Azz,intsteps)+...
                        linemag(k,xstart3,xstop3,ystart3,ystop3,zstart3,zstop3,Axx,Ayy,Azz,intsteps)-...
                        linemag(k,xstart4,xstop4,ystart4,ystop4,zstart4,zstop4,Axx,Ayy,Azz,intsteps)); 

end
Bmag=plotline(B,steps,1);
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
 B(counter,(4:6))= I*(linemag(k,xstart1,xstop1,ystart1,ystop1,zstart1,zstop1,Axx,Ayy,Azz,intsteps)-...
                        linemag(k,xstart2,xstop2,ystart2,ystop2,zstart2,zstop2,Axx,Ayy,Azz,intsteps)+...
                        linemag(k,xstart3,xstop3,ystart3,ystop3,zstart3,zstop3,Axx,Ayy,Azz,intsteps)-...
                        linemag(k,xstart4,xstop4,ystart4,ystop4,zstart4,zstop4,Axx,Ayy,Azz,intsteps)); 

end   

%Bmag=plotline(B,steps,4);
key
grid on
%hold off
%plot(x(:,1),x(:,5),'kh',xy(:,1),xy(:,5),'yo',z(:,1),z(:,5),'c.',y(:,1),y(:,5),'rs')
hold off

toc