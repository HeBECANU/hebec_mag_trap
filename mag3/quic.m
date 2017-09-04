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
intsteps=400;

zah=2*0.5/100;
Rah=0.5/100;
Iah=6;%60; %correct value =118.5 ;146.5 114.5
Nah=30;%23.5

xQ=2*2/100;
RQ=.6/100;%(1.4/2)*(25.4/1000);
IQ=20; %-60
NQ=90;%22.5
%%%% Note 1T = 10^4 G %%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%First CL calc  .%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% over what range? corrdinates in mm
range=15;
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
lstart=5/1000;
lstop=2.00000000001/1000;
%%%%%%%% load parameters
addload =0;
addloadends=0;
ILoad=100;

x1start=-lstart;y1start=lstart;z1start=-30/100;x1stop=-lstop;y1stop=lstop;z1stop=-1/100;
x2start=lstart;y2start=lstart;z2start=-30/100;
x2stop=lstop;y2stop=lstop;z2stop=-1/100;
x3start=-lstart;y3start=-lstart;z3start=-30/100;
x3stop=-lstop;y3stop=-lstop;z3stop=-1/100;
x4start=lstart;y4start=-lstart;z4start=-30/100;
x4stop=lstop;y4stop=-lstop;z4stop=-1/100;

% end coordinates
%x1->x2
x1estart=x1stop;
y1estart=y1stop;
z1estart=z1stop;
x1estop=x2stop;
y1estop=y2stop;
z1estop=z2stop;

%x3->x4
x2estart=x3stop;
y2estart=y3stop;
z2estart=z3stop;
x2estop=x4stop;
y2estop=y4stop;
z2estop=z4stop;

tp=-0.0/100;

% first calculate for z-axis
for ttt=1:steps+1
   
   % assign test points
   Axx=tp;
   Ayy=0;
   Azz=zvec(ttt);
   
 counter=counter+1;
 B(counter,1)=Axx;
 B(counter,2)=Ayy;
 B(counter,3)=Azz-tp;
 B(counter,(4:6))= (Iah*Nah*ahmag(k,Rah,zah,Axx,Ayy,Azz,intsteps))...
    +(IQ*NQ*Qmag(k,RQ,xQ,Axx,Ayy,Azz,intsteps))+[0, 0, 0]; 

if addload ==1
    B(counter,(4:6))=B(counter,(4:6))+(ILoad*Loadmag2(k,x1start,y1start,z1start,x1stop,y1stop,z1stop,...
        x2start,y2start,z2start,x2stop,y2stop,z2stop,x3start,y3start,z3start,x3stop,y3stop,z3stop,...
        x4start,y4start,z4start,x4stop,y4stop,z4stop,Axx,Ayy,Azz,intsteps));
end
if addloadends ==1
    B(counter,(4:6))=B(counter,(4:6))+(ILoad*Loadmagend(k,x1estart,y1estart,z1estart,x1estop,y1estop,z1estop,...
        x2estart,y2estart,z2estart,x2estop,y2estop,z2estop,Axx,Ayy,Azz,intsteps));
end
   

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
 B(counter,(4:6))= (Iah*Nah*ahmag(k,Rah,zah,Axx,Ayy,Azz,intsteps))...
    +(IQ*NQ*Qmag(k,RQ,xQ,Axx,Ayy,Azz,intsteps))+[0, 0, 0]; 

if addload ==1
    B(counter,(4:6))=B(counter,(4:6))+(ILoad*Loadmag2(k,x1start,y1start,z1start,x1stop,y1stop,z1stop,...
        x2start,y2start,z2start,x2stop,y2stop,z2stop,x3start,y3start,z3start,x3stop,y3stop,z3stop,...
        x4start,y4start,z4start,x4stop,y4stop,z4stop,Axx,Ayy,Azz,intsteps));
end
   
if addloadends ==1
    B(counter,(4:6))=B(counter,(4:6))+(ILoad*Loadmagend(k,x1estart,y1estart,z1estart,x1estop,y1estop,z1estop,...
        x2estart,y2estart,z2estart,x2estop,y2estop,z2estop,Axx,Ayy,Azz,intsteps));
end

end
Bmag=plotline(B,steps,2);
counter=0;

% now calculate for x-axis
for ttt=1:steps+1
   
   % assign test points
   Axx=xvec(ttt);
   Ayy=0;
   Azz=0;
   
 counter=counter+1;
 B(counter,1)=Axx;
 B(counter,2)=Ayy;
 B(counter,3)=Azz;
 B(counter,(4:6))= (Iah*Nah*ahmag(k,Rah,zah,Axx,Ayy,Azz,intsteps))...
    +(IQ*NQ*Qmag(k,RQ,xQ,Axx,Ayy,Azz,intsteps))+[0, 0, 0]; 

if addload ==1
    B(counter,(4:6))=B(counter,(4:6))+(ILoad*Loadmag2(k,x1start,y1start,z1start,x1stop,y1stop,z1stop,...
        x2start,y2start,z2start,x2stop,y2stop,z2stop,x3start,y3start,z3start,x3stop,y3stop,z3stop,...
        x4start,y4start,z4start,x4stop,y4stop,z4stop,Axx,Ayy,Azz,intsteps));
end

if addloadends ==1
    B(counter,(4:6))=B(counter,(4:6))+(ILoad*Loadmagend(k,x1estart,y1estart,z1estart,x1estop,y1estop,z1estop,...
        x2estart,y2estart,z2estart,x2estop,y2estop,z2estop,Axx,Ayy,Azz,intsteps));
end   

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
 B(counter,(4:6))= (Iah*Nah*ahmag(k,Rah,zah,Axx,Ayy,Azz,intsteps))...
    +(IQ*NQ*Qmag(k,RQ,xQ,Axx,Ayy,Azz,intsteps))+[0, 0, 0]; 

if addload ==1
    B(counter,(4:6))=B(counter,(4:6))+(ILoad*Loadmag2(k,x1start,y1start,z1start,x1stop,y1stop,z1stop,...
        x2start,y2start,z2start,x2stop,y2stop,z2stop,x3start,y3start,z3start,x3stop,y3stop,z3stop,...
        x4start,y4start,z4start,x4stop,y4stop,z4stop,Axx,Ayy,Azz,intsteps));
end
   
if addloadends ==1
    B(counter,(4:6))=B(counter,(4:6))+(ILoad*Loadmagend(k,x1estart,y1estart,z1estart,x1estop,y1estop,z1estop,...
        x2estart,y2estart,z2estart,x2estop,y2estop,z2estop,Axx,Ayy,Azz,intsteps));
end

end   

%Bmag=plotline(B,steps,4);
key
grid on
%hold off
%plot(x(:,1),x(:,5),'kh',xy(:,1),xy(:,5),'yo',z(:,1),z(:,5),'c.',y(:,1),y(:,5),'rs')
hold off

toc