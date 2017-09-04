%calculates the b-field for one coil


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
Iz=0; %-60
Nz=50;%22.5
zz=2*31/1000;%*(31+(Nah*1.6))/1000;
Rz=0.043;%(1.4/2)*(25.4/1000); Mean radius + a little for heatshrink


%%%%%%%%%%%%%%%%%%%%%%%For coils OD =1.8mm, ID =0.7 mm%%%%%%%%%%%%%
%%%%% R= 21/1000 Ohms/m
Iz2=-10; %-60
Nz2=100;%22.5
zz2=2*40/1000;%*(31+(Nah*1.6))/1000;
Rz2=0.04;%(1.4/2)*(25.4/1000); Mean radius + a little for heatshrink
offset=.05;

%%%%%%%%%%%%%%%%%%%%%%%For coils OD =1.8mm, ID =0.7 mm%%%%%%%%%%%%%
%%%%% R= 21/1000 Ohms/m
Iz3=10; %-60
Nz3=100;%22.5
zz3=2*40/1000;%*(31+(Nah*1.6))/1000;
Rz3=0.04;%(1.4/2)*(25.4/1000); Mean radius + a little for heatshrink
offsety=-.05;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%First CL calc  .%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% over what range? corrdinates in mm
range=20;
xmax=1*range;
xmin=-1*range;
ymax=1*range;
ymin=-1*range;
zmax=1*range;
zmin=-1*range;
steps=10;
[xvec,yvec,zvec]=setupcoords(xmin,xmax,ymin,ymax,zmin,zmax,steps);
counter=0;
B=zeros(21,6);

numm=5;
tpstart=-1/100;
steptp=2*abs(tpstart)/numm;
for sss=1:numm
    tp=tpstart+(steptp*(sss-1))
   
% first calculate for z-axis
for ttt=1:steps+1
   
   % assign test points
   Axx=tp;
   Ayy=0;
   Azz=zvec(ttt);
   
 counter=counter+1;
 B(counter,1)=Axx;
 B(counter,2)=Ayy;
 B(counter,3)=Azz;
 B(counter,(4:6))= (Iz*Nz*onecoilmag(k,Rz,zz,Axx,Ayy,Azz,intsteps))+(Iz2*Nz2*onecoiloffsetx(k,offset,Rz2,zz2,Axx,Ayy,Azz,intsteps))...
     +(Iz3*Nz3*onecoiloffsety(k,offsety,Rz3,zz3,Axx,Ayy,Azz,intsteps));

end
end

quiver3(B(:,1),B(:,2),B(:,3),B(:,4),B(:,5),B(:,6))
pause
subplot(2,1,1)
quiver3(B(:,1),B(:,2),B(:,3),B(:,4),B(:,5),B(:,6))
view(-90,0)

subplot(2,1,2)
quiver3(B(:,1),B(:,2),B(:,3),B(:,4),B(:,5),B(:,6))
view(0,0)

