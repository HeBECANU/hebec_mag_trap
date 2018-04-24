%a program to calculate the field for a 
%QUIC trap. The trap comprises two sets of anti-helmholtz coils
% offset from each other
% Not an extended coil model

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
intsteps=500;
dx=0e-6;
%%%%%%%%%%%%%%%%%First coil set %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
z1=2*.835/100-dx;
R1=.75/100+dx;%4.6228
I1=20; %20
N1=10;%10
%%%%%%%%%%%%%%%%%Second coil set %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
z2=2*0.835/100-dx;
R2=0.7/100+dx;%4.6228
I2=20; %20
N2=20;%26
x2=1.85/100;

%%%% Note 1T = 10^4 G %%%%%%
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
steps=80;

[xvec,yvec,zvec]=setupcoords(xmin,xmax,ymin,ymax,zmin,zmax,steps);
counter=0;
wirecount=0;
B=zeros(21,6);

testpoint=.49/100;

% first calculate for z-axis

for ttt=1:steps+1
% assign test points

    Axx=testpoint;
    Ayy=0;
    Azz=zvec(ttt);

counter=counter+1;
B(counter,1)=Axx;B(counter,2)=Ayy;B(counter,3)=Azz;

B(counter,(4:6))=(I1*N1*ahmag(k,R1,z1,Axx,Ayy,Azz,intsteps))+(I2*N2*ahmagoff(k,x2,R2,z2,Axx,Ayy,Azz,intsteps));
end

Bmag=plotline(B,steps,3);
hold on
counter=0;

% now calculate for y-axis

for ttt=1:steps+1

% assign test points
    Axx=testpoint;
    Ayy=yvec(ttt);
    Azz=0;
counter=counter+1;
B(counter,1)=Axx;B(counter,2)=Ayy;B(counter,3)=Azz;

B(counter,(4:6))=(I1*N1*ahmag(k,R1,z1,Axx,Ayy,Azz,intsteps))+(I2*N2*ahmagoff(k,x2,R2,z2,Axx,Ayy,Azz,intsteps));
end

Bmag=plotline(B,steps,2);
counter=0;

% now calculate for x-axis

for ttt=1:steps+1
% assign test points

    Axx=xvec(ttt)+testpoint;
    Ayy=0;
    Azz=0;
    counter=counter+1;

B(counter,1)=Axx-testpoint;B(counter,2)=Ayy;B(counter,3)=Azz;

B(counter,(4:6))=(I1*N1*ahmag(k,R1,z1,Axx,Ayy,Azz,intsteps))+(I2*N2*ahmagoff(k,x2,R2,z2,Axx,Ayy,Azz,intsteps));

end
Bmag=plotline(B,steps,1);
counter=0;

% now calculate for xy-axis

for ttt=1:steps+1   
% assign test points

Axx=xvec(ttt)+testpoint;
Ayy=yvec(ttt);
Azz=0;
counter=counter+1;

B(counter,1)=Axx-testpoint;B(counter,2)=Ayy;B(counter,3)=Azz;

B(counter,(4:6))=(I1*N1*ahmag(k,R1,z1,Axx,Ayy,Azz,intsteps))+(I2*N2*ahmagoff(k,x2,R2,z2,Axx,Ayy,Azz,intsteps));

end   
Bmag=plotline(B,steps,4);
grid on
hold off

toc