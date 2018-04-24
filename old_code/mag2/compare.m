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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Clover details%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Iclover=150; %130
Rinstart=(1.56)*(25.4/1000);% inner radius in m (distance to the edge of the wire touching the former) 
Routstart=(1.775)*(25.4/1000);%4.775/100;% outer radius in m inner radius in m (distance to the edge of the wire touching the former)
widthofcloverwire=(.04)*(25.4/1000);
numlayerscloverr=4;
numlayerscloverz=4;
totalwiresclover=15;
zcloverstart=(2*(2.9/100))-(numlayerscloverz*widthofcloverwire); % seperation in metres %inner distance (and also the inside of the wire)
clovermat=makewiremat(numlayerscloverr,numlayerscloverz,totalwiresclover,0)
offsetstart=.2585*(25.4/1000);
pause
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Bias details%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Rbiasstart=4.1656/100; %starting radius
Ibias=150; %110
widthofbiaswire=.112/100;
numlayersbiasr=6; %this is the number of coils in the 'r' direction
numlayersbiasz=4; %this is the number of coils in the 'z' direction
totalwiresbias=23;
zbiasstart=(2*2.15/100)-(numlayersbiasz*widthofbiaswire);%inner distance
biasmat=makewiremat(numlayersbiasr,numlayersbiasz,totalwiresbias,0)
pause

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Curvature details%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Rcurvstart=(1.12/2)*(25.4/1000); %starting radius
Icurv=150; %130
widthofcurvwire=(.04)*(25.4/1000);
numlayerscurvr=7;
numlayerscurvz=5;
totalwirescurv=34;
curvmat=makewiremat(numlayerscurvr,numlayerscurvz,totalwirescurv,0)
zcurvstart=(2*2.15/100)-(numlayerscurvz*widthofcurvwire); %inner distance (and also the inside of the wire)
pause

%%%% Note 1T = 10^4 G %%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%First CL calc  .%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% over what range? corrdinates in mm
xmax=10;
xmin=-10;
ymax=10;
ymin=-10;
zmax=10;
zmin=-10;
steps=20;
[xvec,yvec,zvec]=setupcoords(xmin,xmax,ymin,ymax,zmin,zmax,steps);
counter=0;

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
 B(counter,(4:6))=zeros(1,3);
 
 for ss=1:numlayersbiasz
    zbias=zbiasstart+(widthofbiaswire/2)+(widthofbiaswire*(ss-1));
 for tt =1:numlayersbiasr
    Rbias=Rbiasstart+(widthofbiaswire/2)+(widthofbiaswire*(tt-1));
 B(counter,(4:6))=B(counter,(4:6))+(biasmat(tt,ss)*(Ibias*biasmag(k,Rbias,zbias,Axx,Ayy,Azz,intsteps)));
    end
 end

 for ss=1:numlayerscurvz
    zcurv=zcurvstart+(widthofcurvwire/2)+(widthofcurvwire*(ss-1));
 for tt =1:numlayerscurvr
    Rcurv=Rcurvstart+(widthofcurvwire/2)+(widthofcurvwire*(tt-1));
 B(counter,(4:6))=B(counter,(4:6))+(curvmat(tt,ss)*(Icurv*curvmag(k,Rcurv,zcurv,Axx,Ayy,Azz,intsteps))); 
end
 end
 
 for ss=1:numlayerscloverz
    zclover=zcloverstart+(widthofcloverwire/2)+(widthofcloverwire*(ss-1));
 for tt =1:numlayerscloverr
    Rin=Rinstart-(widthofcloverwire/2)-(widthofcloverwire*(tt-1));
    Rout=Routstart+(widthofcloverwire/2)+(widthofcloverwire*(tt-1));
	offset=offsetstart-(widthofcloverwire/2)-(widthofcloverwire*(tt-1));
   B(counter,(4:6))=B(counter,(4:6))+(clovermat(tt,ss)*(Iclover*clovermagcomp(k,Rin,Rout,offset,zclover,Axx,Ayy,Azz,intsteps,0)));
%%%****** Remember this n+1 thing!!!!!!!!!!!!!!!!!!!!!!
end
end

 %+(Icurv*curvmag(k,Rcurv,zcurv,Axx,Ayy,Azz,intsteps)); 
   %-(Iclover*extramag(k,Rin,Rout,z1,Axx,Ayy,Azz,intsteps)); 
%(Iclover*(Nclover+1)*clovermag(k,Rin,Rout,z1,Axx,Ayy,Azz,intsteps,0))...

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
 B(counter,(4:6))=zeros(1,3);
 
for ss=1:numlayersbiasz
    zbias=zbiasstart+(widthofbiaswire/2)+(widthofbiaswire*(ss-1));
 for tt =1:numlayersbiasr
    Rbias=Rbiasstart+(widthofbiaswire/2)+(widthofbiaswire*(tt-1));
 B(counter,(4:6))=B(counter,(4:6))+(biasmat(tt,ss)*(Ibias*biasmag(k,Rbias,zbias,Axx,Ayy,Azz,intsteps)));
    end
 end

 for ss=1:numlayerscurvz
    zcurv=zcurvstart+(widthofcurvwire/2)+(widthofcurvwire*(ss-1));
 for tt =1:numlayerscurvr
    Rcurv=Rcurvstart+(widthofcurvwire/2)+(widthofcurvwire*(tt-1));
 B(counter,(4:6))=B(counter,(4:6))+(curvmat(tt,ss)*(Icurv*curvmag(k,Rcurv,zcurv,Axx,Ayy,Azz,intsteps))); 
end
 end
 
 for ss=1:numlayerscloverz
    zclover=zcloverstart+(widthofcloverwire/2)+(widthofcloverwire*(ss-1));
 for tt =1:numlayerscloverr
    Rin=Rinstart-(widthofcloverwire/2)-(widthofcloverwire*(tt-1));
    Rout=Routstart+(widthofcloverwire/2)+(widthofcloverwire*(tt-1));
	offset=offsetstart-(widthofcloverwire/2)-(widthofcloverwire*(tt-1));
   B(counter,(4:6))=B(counter,(4:6))+(clovermat(tt,ss)*(Iclover*clovermagcomp(k,Rin,Rout,offset,zclover,Axx,Ayy,Azz,intsteps,0)));
%%%****** Remember this n+1 thing!!!!!!!!!!!!!!!!!!!!!!
end
end
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
 B(counter,(4:6))=zeros(1,3);
 
 for ss=1:numlayersbiasz
    zbias=zbiasstart+(widthofbiaswire/2)+(widthofbiaswire*(ss-1));
 for tt =1:numlayersbiasr
    Rbias=Rbiasstart+(widthofbiaswire/2)+(widthofbiaswire*(tt-1));
 B(counter,(4:6))=B(counter,(4:6))+(biasmat(tt,ss)*(Ibias*biasmag(k,Rbias,zbias,Axx,Ayy,Azz,intsteps)));
    end
 end

 for ss=1:numlayerscurvz
    zcurv=zcurvstart+(widthofcurvwire/2)+(widthofcurvwire*(ss-1));
 for tt =1:numlayerscurvr
    Rcurv=Rcurvstart+(widthofcurvwire/2)+(widthofcurvwire*(tt-1));
 B(counter,(4:6))=B(counter,(4:6))+(curvmat(tt,ss)*(Icurv*curvmag(k,Rcurv,zcurv,Axx,Ayy,Azz,intsteps))); 
end
 end
 
 for ss=1:numlayerscloverz
    zclover=zcloverstart+(widthofcloverwire/2)+(widthofcloverwire*(ss-1));
 for tt =1:numlayerscloverr
    Rin=Rinstart-(widthofcloverwire/2)-(widthofcloverwire*(tt-1));
    Rout=Routstart+(widthofcloverwire/2)+(widthofcloverwire*(tt-1));
	offset=offsetstart-(widthofcloverwire/2)-(widthofcloverwire*(tt-1));
   B(counter,(4:6))=B(counter,(4:6))+(clovermat(tt,ss)*(Iclover*clovermagcomp(k,Rin,Rout,offset,zclover,Axx,Ayy,Azz,intsteps,0)));
%%%****** Remember this n+1 thing!!!!!!!!!!!!!!!!!!!!!!
end
end
end

Bmag=plotline(B,steps,4);
key
hold off
toc