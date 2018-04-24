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
k=(Mu0)/(4*pi);
intsteps=200;

dz = -.0;
addedbias=0;%-15e-4;
gasketthickness=.05;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Clover details%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Rclover=2*.135;
Vclover=50;
Iclover=150; %150
Rinstart=(.780)*(25.4/1000);% inner radius in m (distance to the edge of the wire touching the former) 
Rinstop=(.557)*(25.4/1000);
Routstart=(1.775)*(25.4/1000);%4.775/100;% outer radius in m inner radius in m (distance to the edge of the wire touching the former)
Routstop=(2.0)*(25.4/1000);
widthofcloverwire=.112/100;%(.04)*(25.4/1000);
numlayerscloverr=3;%3
numlayerscloverz=5;%4
totalwiresclover=14; %mean n 10
zcloverstart=(2*((2.66251+gasketthickness+dz)/100)); % seperation in metres %inner distance (and also the inside of the wire)
zcloverstop=(2*((3.33+gasketthickness+dz)/100));%3.13749
clovermat=makewiremat(numlayerscloverr,numlayerscloverz,totalwiresclover,0)
offsetstart=.2585*(25.4/1000);
offsetstop=(.067/2)*(25.4/1000);


pause
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Bias details%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Rbiasandcurv=2*(.18+.09);
Ibias=115;
Rbiasstart=(3.104/2)*(25.4/1000); %starting radius
Rbiasstop=(3.669/2)*(25.4/1000); % stopping radius
widthofbiaswire=.112/100;%(.04)*(25.4/1000);
numlayersbiasr=6; %this is the number of coils in the 'r' direction
numlayersbiasz=5; %this is the number of coils in the 'z' direction
totalwiresbias=30;
zbiasstart=(2*(1.788+dz)/100);%inner distance
zbiasstop=(2*(2.512+dz)/100);%outer distance
biasmat=makewiremat(numlayersbiasr,numlayersbiasz,totalwiresbias,0)
biastheta=0; % extra bit of turn to get wire out
pause

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Curvature details%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Rcurvstart=(1.12/2)*(25.4/1000); %starting radius
Rcurvstop=(1.765/2)*(25.4/1000) %starting radius
Icurv=115; %157.1, 132.7 freq@132.7 = 106 Hz, 127.6
widthofcurvwire=.112/100;%(.04)*(25.4/1000);
numlayerscurvr=7;%7
numlayerscurvz=5;%5
totalwirescurv=34;
curvmat=makewiremat(numlayerscurvr,numlayerscurvz,totalwirescurv,0)
zcurvstart=(2*(1.788+dz)/100); %inner distance (and also the inside of the wire)
zcurvstop=(2*(2.512+dz)/100);
curvtheta=0;
pause

%%%% Note 1T = 10^4 G %%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%First CL calc  .%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% over what range? corrdinates in mm
xmax=15;
xmin=-15;
ymax=15;
ymin=-15;
zmax=15;
zmin=-15;
steps=40;
[xvec,yvec,zvec]=setupcoords(xmin,xmax,ymin,ymax,zmin,zmax,steps);
counter=0;

 %allow space for each wire 
biasgapR=findgap(numlayersbiasr,widthofbiaswire,Rbiasstop,Rbiasstart); 
biasgapz=findgap(numlayersbiasz,widthofbiaswire,zbiasstop,zbiasstart); 
curvgapR=findgap(numlayerscurvr,widthofcurvwire,Rcurvstop,Rcurvstart); 
curvgapz=findgap(numlayerscurvz,widthofcurvwire,zcurvstop,zcurvstart); 
clovgapRin=findgap(numlayerscloverr,widthofcloverwire,Rinstart,Rinstop); 
clovgapRout=findgap(numlayerscloverr,widthofcloverwire,Routstop,Routstart); 
clovgapoffset=findgap(numlayerscloverr,widthofcloverwire,offsetstart,offsetstop); 
clovergapz=findgap(numlayerscloverz,widthofcloverwire,zcloverstop,zcloverstart); 

 %ammend z,R starts to allow for initial gap
 zbiasstart=zbiasstart+biasgapz;
 Rbiasstart=Rbiasstart+biasgapR;
 zcurvstart=zcurvstart+curvgapz;
 Rcurvstart=Rcurvstart+curvgapR;
 zcloverstart =zcloverstart+clovergapz;
 Rinstart=Rinstart+clovgapRin;
 Routstart=Routstart+clovgapRout;
 offsetstart=offsetstart+clovgapoffset;

% calculate for zx-axis
for ttt=1:steps+1
   counter=0;
   for sss=1:steps+1

   % assign test points
   Axx=xvec(sss);
   Ayy=0;
   Azz=zvec(ttt);
   
 counter=counter+1;
 B(counter,1)=Axx;
 B(counter,2)=Ayy;
 B(counter,3)=Azz;
 B(counter,(4:6))=zeros(1,3);
  
 for ss=1:numlayersbiasz
    zbias=zbiasstart+(widthofbiaswire/2)+(widthofbiaswire*(ss-1))+(2*biasgapz*(ss-1));
 for tt =1:numlayersbiasr
    Rbias=Rbiasstart+(widthofbiaswire/2)+(widthofbiaswire*(tt-1))+(2*biasgapR*(tt-1));
 B(counter,(4:6))=B(counter,(4:6))+(biasmat(tt,ss)*(Ibias*biasmag(k,Rbias,zbias,Axx,Ayy,Azz,intsteps)));
    end
 end
 

 if biastheta ~=0
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %%Add in extra half wrap for bias coil to get out%%%%%%%%%%%%%%%%%%
 B(counter,(4:6))=B(counter,(4:6))+(Ibias*biasmagextra(k,Rbias,zbias,Axx,Ayy,Azz,intsteps,biastheta));
end

 for ss=1:numlayerscurvz
    zcurv=zcurvstart+(widthofcurvwire/2)+(widthofcurvwire*(ss-1))+(2*curvgapz*(ss-1));
 for tt =1:numlayerscurvr
    Rcurv=Rcurvstart+(widthofcurvwire/2)+(widthofcurvwire*(tt-1))+(2*curvgapR*(tt-1));
 B(counter,(4:6))=B(counter,(4:6))+(curvmat(tt,ss)*(Icurv*curvmag(k,Rcurv,zcurv,Axx,Ayy,Azz,intsteps))); 
end
end
 %B(counter,(4:6))=B(counter,(4:6))+(2*(Icurv*andymag(k,Rcurv,zcurv,Axx,Ayy,Azz,intsteps))); 

if curvtheta~= 0
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %%Add in extra half wrap for curv coil to get out%%%%%%%%%%%%%%%%%%
 B(counter,(4:6))=B(counter,(4:6))+(Icurv*curvmagextra(k,Rcurv,zcurv,Axx,Ayy,Azz,intsteps,curvtheta));
end

for ss=1:numlayerscloverz
    zclover=zcloverstart+(widthofcloverwire/2)+(widthofcloverwire*(ss-1))+(2*clovergapz*(ss-1));
 for tt =1:numlayerscloverr
    Rin=Rinstart-(widthofcloverwire/2)-(widthofcloverwire*(tt-1))-(2*clovgapRin*(tt-1));
    Rout=Routstart+(widthofcloverwire/2)+(widthofcloverwire*(tt-1))+(2*clovgapRout*(tt-1));
    offset=offsetstart-(widthofcloverwire/2)-(widthofcloverwire*(tt-1))-(2*clovgapoffset*(tt-1));
   B(counter,(4:6))=B(counter,(4:6))+(clovermat(tt,ss)*(Iclover*clovermagcomp(k,Rin,Rout,offset,zclover,Axx,Ayy,Azz,intsteps,0,totalwiresclover,(tt*ss))));
end
end
   B(counter,(4:6))=B(counter,(4:6))+(Iclover*clovermagcompextra(k,Rin,Rout,offset,zclover,Axx,Ayy,Azz,intsteps,0)); % this adds the extra outer arc
end
A= sqrt(sum(B(:,(4:6)).^2,2))/1e-4;
Bmag(:,ttt)=A;
end


surf(xvec,zvec,Bmag)