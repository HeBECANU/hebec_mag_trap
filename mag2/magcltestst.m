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

zcloverstart=2*2.9/100; % seperation in metres %inner distance (and also the inside of the wire)
Iclover=130; %130
Rinstart=1.69799/100;% inner radius in m (distance to the edge of the wire touching the former) 
Routstart=4.79425/100;%4.775/100;% outer radius in m inner radius in m (distance to the edge of the wire touching the former)
widthofcloverwire=1/100000000;
numwirescloverr=3;
numwirescloverz=2;
totalwiresclover=numwirescloverr*numwirescloverz;
offsetstart=1/100000000;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Bias details%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

zbiasstart=2*2.15/100;%inner distance
Rbiasstart=4.301/100; %starting radius
Ibias=110; %110
widthofbiaswire=1/100000000;
numwiresbiasr=2;
numwiresbiasz=3;
totalwiresbias=numwiresbiasr*numwiresbiasz;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Curvature details%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

zcurvstart=2*2.15/100; %inner distance (and also the inside of the wire)
Rcurvstart=1.8319/100; %starting radius
Icurv=130; %130
widthofcurvwire=1/10000000;
numwirescurvr=3;
numwirescurvz=2;
totalwirescurv=numwirescurvr*numwirescurvz;

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
 
 for ss=1:numwiresbiasz
    zbias=zbiasstart+(widthofbiaswire/2)+(widthofbiaswire*ss);
 for tt =1:numwiresbiasr
    Rbias=Rbiasstart+(widthofbiaswire/2)+(widthofbiaswire*tt);
 B(counter,(4:6))=B(counter,(4:6))+(Ibias*biasmag(k,Rbias,zbias,Axx,Ayy,Azz,intsteps));
    end
 end

 for ss=1:numwirescurvz
    zcurv=zcurvstart+(widthofcurvwire/2)+(widthofcurvwire*ss);
 for tt =1:numwirescurvr
    Rcurv=Rcurvstart+(widthofcurvwire/2)+(widthofcurvwire*tt);
 B(counter,(4:6))=B(counter,(4:6))+(Icurv*curvmag(k,Rcurv,zcurv,Axx,Ayy,Azz,intsteps)); 
end
 end
 
 for ss=1:numwirescloverz
    zclover=zcloverstart+(widthofcloverwire/2)+(widthofcloverwire*ss);
 for tt =1:numwirescloverr
    Rin=Rinstart-(widthofcloverwire/2)-(widthofcloverwire*tt);
    Rout=Routstart+(widthofcloverwire/2)+(widthofcloverwire*tt);
	offset=offsetstart-(widthofcloverwire/2)-(widthofcloverwire*tt);
   B(counter,(4:6))=B(counter,(4:6))+(Iclover*clovermagcomp(k,Rin,Rout,offset,zclover,Axx,Ayy,Azz,intsteps,0));
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
 
 for ss=1:numwiresbiasz
    zbias=zbiasstart+(widthofbiaswire/2)+(widthofbiaswire*ss);
 for tt =1:numwiresbiasr
    Rbias=Rbiasstart+(widthofbiaswire/2)+(widthofbiaswire*tt);
 B(counter,(4:6))=B(counter,(4:6))+(Ibias*biasmag(k,Rbias,zbias,Axx,Ayy,Azz,intsteps));
    end
 end

 for ss=1:numwirescurvz
    zcurv=zcurvstart+(widthofcurvwire/2)+(widthofcurvwire*ss);
 for tt =1:numwirescurvr
    Rcurv=Rcurvstart+(widthofcurvwire/2)+(widthofcurvwire*tt);
 B(counter,(4:6))=B(counter,(4:6))+(Icurv*curvmag(k,Rcurv,zcurv,Axx,Ayy,Azz,intsteps)); 
end
 end
 
 for ss=1:numwirescloverz
    zclover=zcloverstart+(widthofcloverwire/2)+(widthofcloverwire*ss);
 for tt =1:numwirescloverr
    Rin=Rinstart-(widthofcloverwire/2)-(widthofcloverwire*tt);
    Rout=Routstart+(widthofcloverwire/2)+(widthofcloverwire*tt);
	offset=offsetstart-(widthofcloverwire/2)-(widthofcloverwire*tt);
B(counter,(4:6))=B(counter,(4:6))+(Iclover*clovermagcomp(k,Rin,Rout,offset,zclover,Axx,Ayy,Azz,intsteps,0));
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
 
 for ss=1:numwiresbiasz
    zbias=zbiasstart+(widthofbiaswire/2)+(widthofbiaswire*ss);
 for tt =1:numwiresbiasr
    Rbias=Rbiasstart+(widthofbiaswire/2)+(widthofbiaswire*tt);
 B(counter,(4:6))=B(counter,(4:6))+(Ibias*biasmag(k,Rbias,zbias,Axx,Ayy,Azz,intsteps));
    end
 end

 for ss=1:numwirescurvz
    zcurv=zcurvstart+(widthofcurvwire/2)+(widthofcurvwire*ss);
 for tt =1:numwirescurvr
    Rcurv=Rcurvstart+(widthofcurvwire/2)+(widthofcurvwire*tt);
 B(counter,(4:6))=B(counter,(4:6))+(Icurv*curvmag(k,Rcurv,zcurv,Axx,Ayy,Azz,intsteps)); 
end
 end
 
 for ss=1:numwirescloverz
    zclover=zcloverstart+(widthofcloverwire/2)+(widthofcloverwire*ss);
 for tt =1:numwirescloverr
    Rin=Rinstart-(widthofcloverwire/2)-(widthofcloverwire*tt);
    Rout=Routstart+(widthofcloverwire/2)+(widthofcloverwire*tt);
	offset=offsetstart-(widthofcloverwire/2)-(widthofcloverwire*tt);
B(counter,(4:6))=B(counter,(4:6))+(Iclover*clovermagcomp(k,Rin,Rout,offset,zclover,Axx,Ayy,Azz,intsteps,0));
%%%****** Remember this n+1 thing!!!!!!!!!!!!!!!!!!!!!!
end
end

end

Bmag=plotline(B,steps,4);

hold off
toc