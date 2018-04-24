%calculates the magnetic field of an extended double wound solenoid
% as in Opat paper

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

%%%%%First solenoid parameters
I=3; %-60
N=100;%22.5
R=.3/100;
T1=.5/1000; %Wire thickness
gap=.9/1000; %gap from edge of T1 to edge of T2, i.e. centre to centre =T1/2+T2/2+gap!!!!
H=2/1000;
%%%%%second solenoid parameters
I2=3;
T2=.5/1000; %Wire thickness
zstart=-(N*(T1))/2;%start of solenoid 
%%%%Note 1T = 10^4 G %%%%%%

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
steps=50;
[xvec,yvec,zvec]=setupcoords(xmin,xmax,ymin,ymax,zmin,zmax,steps);
counter=0;
B=zeros(21,6);
ab=0e-4;
tp=-0/100;

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
 %%%%%Extended integration %%%%%%
for tt =1:N
    z1temp=zstart+(T1/2)+((T1+T2+(2*gap))*(tt-1));
    B(counter,(4:6))=B(counter,(4:6))+(I*coilmag2(k,R+T1/2,z1temp,Axx,Ayy,Azz,intsteps));  
 %   B(counter,(4:6))=B(counter,(4:6))+(I*linemag(k,1e-12,R+(T1/2)+H,z1temp,0,R+(T1/2)+H,z1temp+gap,Axx,Ayy,Azz,intsteps));
end
for tt =1:N
    z1temp=zstart+(T2/2+T1+gap)+((T2+T1+(2*gap))*(tt-1));
    B(counter,(4:6))=B(counter,(4:6))-(I2*coilmag2(k,R+(T2/2),z1temp,Axx,Ayy,Azz,intsteps));  
  %   B(counter,(4:6))=B(counter,(4:6))-(I2*linemag(k,1e-12,R+(T2/2)+H,z1temp,0,R+(T2/2)+H,z1temp+gap,Axx,Ayy,Azz,intsteps));
end
B(counter,(4:6))=B(counter,(4:6))+[ab 0 0];
end
Bmag=plotline(B,steps,3);
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
 B(counter,(4:6))=zeros(1,3);
 %%%%%Extended integration %%%%%%
 %%%%%Extended integration %%%%%%
for tt =1:N
    z1temp=zstart+(T1/2)+((T1+T2+(2*gap))*(tt-1));
    B(counter,(4:6))=B(counter,(4:6))+(I*coilmag2(k,R+T1/2,z1temp,Axx,Ayy,Azz,intsteps));  
 %   B(counter,(4:6))=B(counter,(4:6))+(I*linemag(k,1e-12,R+(T1/2)+H,z1temp,0,R+(T1/2)+H,z1temp+gap,Axx,Ayy,Azz,intsteps));
end
for tt =1:N
     z1temp=zstart+(T2/2+T1+gap)+((T2+T1+(2*gap))*(tt-1));
    B(counter,(4:6))=B(counter,(4:6))-(I2*coilmag2(k,R+(T2/2),z1temp,Axx,Ayy,Azz,intsteps));  
  %   B(counter,(4:6))=B(counter,(4:6))-(I2*linemag(k,1e-12,R+(T2/2)+H,z1temp,0,R+(T2/2)+H,z1temp+gap,Axx,Ayy,Azz,intsteps));
end
B(counter,(4:6))=B(counter,(4:6))+[ab 0 0];
end
Bmag=plotline(B,steps,2);
counter=0;

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
 B(counter,(4:6))=zeros(1,3);
 %%%%%Extended integration %%%%%%
for tt =1:N
    z1temp=zstart+(T1/2)+((T1+T2+(2*gap))*(tt-1));
    B(counter,(4:6))=B(counter,(4:6))+(I*coilmag2(k,R+T1/2,z1temp,Axx,Ayy,Azz,intsteps));  
  %  B(counter,(4:6))=B(counter,(4:6))+(I*linemag(k,1e-12,R+(T1/2)+H,z1temp,0,R+(T1/2)+H,z1temp+gap,Axx,Ayy,Azz,intsteps));
end
for tt =1:N
    z1temp=zstart+(T2/2+T1+gap)+((T2+T1+(2*gap))*(tt-1));
    B(counter,(4:6))=B(counter,(4:6))-(I2*coilmag2(k,R+(T2/2),z1temp,Axx,Ayy,Azz,intsteps));  
   %  B(counter,(4:6))=B(counter,(4:6))-(I2*linemag(k,1e-12,R+(T2/2)+H,z1temp,0,R+(T2/2)+H,z1temp+gap,Axx,Ayy,Azz,intsteps));
end
B(counter,(4:6))=B(counter,(4:6))+[ab 0 0];
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
 B(counter,(4:6))=zeros(1,3);
 %%%%%Extended integration %%%%%%
for tt =1:N
    z1temp=zstart+(T1/2)+((T1+T2+(2*gap))*(tt-1));
    B(counter,(4:6))=B(counter,(4:6))+(I*coilmag2(k,R+T1/2,z1temp,Axx,Ayy,Azz,intsteps));  
  %  B(counter,(4:6))=B(counter,(4:6))+(I*linemag(k,0,R+(T1/2)+H,z1temp,0,R+(T1/2)+H,z1temp+gap,Axx,Ayy,Azz,intsteps));
end
for tt =1:N
     z1temp=zstart+(T2/2+T1+gap)+((T2+T1+(2*gap))*(tt-1));
    B(counter,(4:6))=B(counter,(4:6))-(I2*coilmag2(k,R+(T2/2),z1temp,Axx,Ayy,Azz,intsteps));  
 %    B(counter,(4:6))=B(counter,(4:6))-(I2*linemag(k,0,R+(T2/2)+H,z1temp,0,R+(T2/2)+H,z1temp+gap,Axx,Ayy,Azz,intsteps));
end

B(counter,(4:6))=B(counter,(4:6))+[ab 0 0];
end   

%Bmag=plotline(B,steps,4);
key
grid on
%hold off
%plot(x(:,1),x(:,5),'kh',xy(:,1),xy(:,5),'yo',z(:,1),z(:,5),'c.',y(:,1),y(:,5),'rs')
hold off

toc