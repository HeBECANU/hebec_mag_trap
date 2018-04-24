%a program to calculate the field for a 
%QUIC trap. The trap comprises two sets of anti-helmholtz coils
% offset from each other
% This is an extended coil model
% this program rectifies the 'z' extended error i.e. factor of two!!!!!
% also includes the extra 1/3 turns
% Same as antihelmquicmodified2 accept it uses helixes instead of loops
%%also includes all the input and output wires

%same as antihelmquichelix accept here we use rectangular wire and
%integrate in the long (radial) direction of the wire

%same as antihelmquichelixtest accept here each axis integration is in a
%loop so that the code is considerably simpler.

%same as ahqh, used for testing purposes.
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

alpha=17e-6; % linear thermal expansion coefficient/ degree, for copper
dT=0;%ten degree change!!
dx=alpha*dT*(0.5/1000);
%%%%%%%%%%%%%%%%%First coil set %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%z-dist. = 4.5mm (dist from centre to glass face)+ 2.6 (glass thickness)+
%           0.8 (thickness of metal former) +.25( thickness of rim)
%           +0.2(gap between glass and former)
zdd=7.5+2.7+0.8+.25+.2;%5+2.6+0.75+0.25+.2;
z1=2*zdd/1000; % distance just before the first wire i.e. ->|o
R1=1.0/100;%4.6228
I1=40*(1-0); %20 -> tight, 15 -> high bias for experiment
N1=14;%23.5
thick1z=0.56/1000+0/1000; % wire thickness in the z direction
thick1r=0.56/1000 %% wire thickness in the r direction
stepsr =1; % integration steps in the r-direction
dr1=thick1r/stepsr; %infintessimal in the r dir.
numz1s=14;

zdd=7.5+2.7+0.8+.25+.2;%5+2.6+0.75+0.25+.2;
%%%%%%%%%%%%%%%%%Second coil set %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
z2=2*zdd/1000;
R2=.7/100;%4.6228
I2=-0*(1+0); %20 -> tight, 30 -> high bias for experiment
N2=10;%23.5
x2=1.85/100; % use 1.95 for higher trap!
thick2z=0.56/1000+0/1000; % wire thickness in the z direction
thick2r=0.56/1000 %% wire thickness in the r direction
stepsr =1; % integration steps in the r-direction
dr2=thick2r/stepsr; %infintessimal in the r dir.
numz2s=10;

%%%%%%%%%%%Line set%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Coordinates of wires%%%
z3=2*zdd/1000;
xs=1.3/100;ys=-1.6/100;xf=xs;yf=-1*ys;zoff=3/100;

N3=18;
I3=-40;
thick3z=.56/1000+0/1000; % wire thickness in the z direction
thick3r=.56/1000 %% wire thickness in the r direction
stepsr =1; % integration steps in the r-direction
dr3=thick3r/stepsr; %infintessimal in the r dir.
numz3s=18;

%%%%%%%%%%%second Line set%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Coordinates of wires%%%
z3=2*zdd/1000;
xs4=1.55/100;ys4=-1.6/100;xf4=xs4;yf4=-1*ys4;zoff=3/100;

N4=18;
I4=-40;
thick3z=.56/1000+0/1000; % wire thickness in the z direction
thick3r=.56/1000 %% wire thickness in the r direction
stepsr =1; % integration steps in the r-direction
dr3=thick3r/stepsr; %infintessimal in the r dir.
numz4s=18;
%assume the wire coating is around 50 microns thick

%%%%%%%%%%%%%%%%%%%%%%%Work out power consumed


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
addwires=0;

testpointx=0.3/100;%0.5/100 -> tight 0.23/100 -> high bias experiment
testpointy=-0.0/100;
addedbias=0e-4;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%What axis do we want plotted%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       [x y z xy]
axiss = [1 0 0 0]

for ll=1:length(axiss)
for ttt=1:steps+1
% check if this axis is wanted %
    if axiss(ll)==0
        break
    end
    wirecount1=0;
    wirecount2=0;
    wirecount3=0;
    wirecount4=0;
% assign test points
if ll == 1
    Axx=xvec(ttt)+testpointx; Ayy=testpointy; Azz=0;
end
if ll == 2
    Axx=testpointx; Ayy=yvec(ttt)+testpointy; Azz=0;
end
if ll == 3
    Axx=testpointx; Ayy=testpointy; Azz=zvec(ttt);
end
if ll == 4
    Axx=xvec(ttt)+testpointx; Ayy=yvec(ttt)+testpointy; Azz=0;
end
    counter=counter+1;

B(counter,1)=Axx;
B(counter,2)=Ayy;
B(counter,3)=Azz;
B(counter,(4:6))=zeros(1,3);

if ll == 1 || ll == 4
   B(counter,1)=Axx-testpointx;
    B(counter,2)=Ayy-testpointy;
    B(counter,3)=Azz;
    B(counter,(4:6))=zeros(1,3); 
end
if I1~= 0
%%%%%%%%Do wires first%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
beta=atan(1.1/.5);
wxstart = (R1+(1.5*thick1r))*cos(beta); wxstop=wxstart;
wystart = (R1+(1.5*thick1r))*sin(beta); wystop=wystart;
wzstart = z1/2+12.5/1000; wzstop = z1/2 +thick1z/2;
%%starting with Quad wires%%%%%%%%%%%%%%%%
%%quad wire in starting from pos x,y and z
if addwires ==1
B(counter,(4:6))=B(counter,(4:6))+(I1*linemag(k,wxstart,wxstop,wystart,wystop,wzstart,wzstop,Axx,Ayy,Azz,intsteps)); 
%%quad wire out starting from neg x,y and pos z
B(counter,(4:6))=B(counter,(4:6))+(I1*linemag(k,wxstart,wxstop,-wystart,-wystop,wzstop,wzstart,Axx,Ayy,Azz,intsteps)); 
%%quad wire in starting from neg x,y and neg z
B(counter,(4:6))=B(counter,(4:6))+(I1*linemag(k,wxstart,wxstop,-wystart,-wystop,-wzstart,-wzstop,Axx,Ayy,Azz,intsteps)); 
%%quad wire out starting from pos x,y and neg z
B(counter,(4:6))=B(counter,(4:6))+(I1*linemag(k,wxstart,wxstop,wystart,wystop,-wzstop,-wzstart,Axx,Ayy,Azz,intsteps)); 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%Extended integration for on-axis coil%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ss=1:numz1s
    if (wirecount1==N1)
    break
    end
    z1temp=z1/2+(thick1z/2)+(thick1z*(ss-1)); % distance to the centre of each wire 
for tt =1:stepsr
    R1temp=R1+(dr1/2)+(dr1*(tt-1));
%%break this into two lots one from theta =beta to theta =0 and then one
%%from 2*pi to beta %%%%%%%%%%%%%%%%%%%%%%%%%
   B(counter,(4:6))=B(counter,(4:6))+((I1/stepsr)*ahmaghelixextra(k,R1temp,beta,0,2*z1temp,thick1z,Axx,Ayy,Azz,intsteps));
   B(counter,(4:6))=B(counter,(4:6))+((I1/stepsr)*ahmaghelixextra(k,R1temp,2*pi,beta,(2*z1temp)+(beta/(pi)*thick1z),thick1z,Axx,Ayy,Azz,intsteps));
    wirecount1=wirecount1+1;
        if (wirecount1==N1)
        break
        end
end
end

%%%%%%%%%Add in the extra bit for the wire to make it out%%%%%%%%%%%%%%%%%%%%%%
%%%do it in two bits, so we don't have theta changing sign%%%%%%%%
for tt =1:stepsr
    R1temp=R1+(dr1/2)+(dr1*(tt-1));
   B(counter,(4:6))=B(counter,(4:6))+((I1/stepsr)*ahmaghelixextra(k,R1temp,beta,0,2*(z1temp+thick1z),thick1z,Axx,Ayy,Azz,intsteps));
   B(counter,(4:6))=B(counter,(4:6))+((I1/stepsr)*ahmaghelixextra(k,R1temp,2*pi,((2*pi)-beta),(2*(z1temp+thick1z))+(beta/(pi)*thick1z),thick1z,Axx,Ayy,Azz,intsteps));
end
end
if I2~=0
%%%%Input Ioffe wires%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
beta=atan(1.1/.5);
wxstart = x2-((R2+(1.5*thick2z))*cos(beta)); wxstop=wxstart;
wystart = (R2+(1.5*thick2z))*sin(beta); wystop=wystart;
wzstart = z1/2+12.5/1000; wzstop = z1/2 +thick2z/2;
%%Ioffe wire in starting from pos x,y and z
if addwires ==1
B(counter,(4:6))=B(counter,(4:6))+(I2*linemag(k,wxstart,wxstop,wystart,wystop,wzstart,wzstop,Axx,Ayy,Azz,intsteps)); 
%%quad wire out starting from neg x,y and pos z
B(counter,(4:6))=B(counter,(4:6))+(I2*linemag(k,wxstart,wxstop,-wystart,-wystop,wzstop,wzstart,Axx,Ayy,Azz,intsteps)); 
%%quad wire in starting from neg x,y and neg z
B(counter,(4:6))=B(counter,(4:6))+(I2*linemag(k,wxstart,wxstop,-wystart,-wystop,-wzstart,-wzstop,Axx,Ayy,Azz,intsteps)); 
%%quad wire out starting from pos x,y and neg z
B(counter,(4:6))=B(counter,(4:6))+(I2*linemag(k,wxstart,wxstop,wystart,wystop,-wzstop,-wzstart,Axx,Ayy,Azz,intsteps)); 
end
   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%%%%Extended integration for off-axis coil%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ss=1:numz2s
    if (wirecount2==N2)
    break
    end
    z2temp=z2/2+(thick2z/2)+(thick2z*(ss-1));   
    
    for tt =1:stepsr
    R2temp=R2+(dr2/2)+(dr2*(tt-1));  
   B(counter,(4:6))=B(counter,(4:6))+((I2/stepsr)*ahmagoffhelix(k,x2,R2temp,pi-beta,2*pi,2*z2temp,thick2z,Axx,Ayy,Azz,intsteps));
   B(counter,(4:6))=B(counter,(4:6))+((I2/stepsr)*ahmagoffhelix(k,x2,R2temp,0,pi-beta,(2*z2temp)+((pi+beta)/pi*thick2z),thick2z,Axx,Ayy,Azz,intsteps));
   wirecount2=wirecount2+1;
        if (wirecount2==N2)
        break
        end
    end
end  
%%%%%%%%%Add in the extra bit for the wire to make it%%%%%%%%%%%%%%%%%%%%%%
for tt =1:stepsr
    R2temp=R2+(dr2/2)+(dr2*(tt-1));
 B(counter,(4:6))=B(counter,(4:6))+((I2/stepsr)*ahmagoffhelix(k,x2,R2temp,pi-beta,pi+beta,2*(z2temp+thick2z),thick2z,Axx,Ayy,Azz,intsteps));
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%Compare line configuration%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if I3~=0
for ss=1:numz3s
    if (wirecount3==N3)
    break
    end
    z3temp=z3/2+(thick3z/2)+(thick3z*(ss-1)); 
    yftemp=yf-(thick3z*(ss-1)); 
    ystemp=ys+(thick3z*(ss-1)); 
    zotemp=zoff-(thick3z*(ss-1)); 
    %%coils on one 'z' side%
   B(counter,(4:6))=B(counter,(4:6))+(I3*linemag(k,xs,xf,yftemp,ystemp,z3temp,z3temp,Axx,Ayy,Azz,intsteps)); %bit that does correct field
 B(counter,(4:6))=B(counter,(4:6))+(I3*linemag(k,xf,xf,ystemp,ystemp,z3temp,z3temp+zotemp,Axx,Ayy,Azz,intsteps)); % bit that heads out
  B(counter,(4:6))=B(counter,(4:6))+(I3*linemag(k,xs,xs,yftemp,yftemp,z3temp+zotemp,z3temp,Axx,Ayy,Azz,intsteps)); % bit that heads in
  B(counter,(4:6))=B(counter,(4:6))+(I3*linemag(k,xs,xf,ystemp,yftemp,z3temp+zotemp,z3temp+zotemp,Axx,Ayy,Azz,intsteps)); % bit does incorrect field
      
   %%coils on one '-z' side%
   B(counter,(4:6))=B(counter,(4:6))+(I3*linemag(k,xs,xf,ystemp,yftemp,-z3temp,-z3temp,Axx,Ayy,Azz,intsteps)); %bit that does correct field
  B(counter,(4:6))=B(counter,(4:6))+(I3*linemag(k,xf,xf,yftemp,yftemp,-z3temp,-(z3temp+zotemp),Axx,Ayy,Azz,intsteps)); % bit that heads out
  B(counter,(4:6))=B(counter,(4:6))+(I3*linemag(k,xs,xf,ystemp,ystemp,-(z3temp+zotemp),-z3temp,Axx,Ayy,Azz,intsteps)); % bit that heads in
  B(counter,(4:6))=B(counter,(4:6))+(I3*linemag(k,xs,xf,yftemp,ystemp,-(z3temp+zotemp),-(z3temp+zotemp),Axx,Ayy,Azz,intsteps)); % bit does incorrect field
      
   wirecount3=wirecount3+1;
        if (wirecount3==N3)
        break
        end

end  
end  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%second line configuration%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if I4~=0
for ss=1:numz4s
    if (wirecount4==N4)
    break
    end
    z4temp=z3/2+(thick3z/2)+(thick3z*(ss-1)); 
    yftemp=yf4-(thick3z*(ss-1)); 
    ystemp=ys4+(thick3z*(ss-1)); 
    %%coils on one 'z' side%
   B(counter,(4:6))=B(counter,(4:6))+(I4*linemag(k,xs4,xf4,yftemp,ystemp,z4temp,z4temp,Axx,Ayy,Azz,intsteps)); %bit that does correct field
 B(counter,(4:6))=B(counter,(4:6))+(I4*linemag(k,xf4,xf4,ystemp,ystemp,z4temp,z4temp+zoff,Axx,Ayy,Azz,intsteps)); % bit that heads out
  B(counter,(4:6))=B(counter,(4:6))+(I4*linemag(k,xs4,xs4,yftemp,yftemp,z4temp+zoff,z4temp,Axx,Ayy,Azz,intsteps)); % bit that heads in
  B(counter,(4:6))=B(counter,(4:6))+(I4*linemag(k,xs4,xf4,ystemp,yftemp,z4temp+zoff,z4temp+zoff,Axx,Ayy,Azz,intsteps)); % bit does incorrect field
      
   %%coils on one '-z' side%
   B(counter,(4:6))=B(counter,(4:6))+(I4*linemag(k,xs4,xf4,ystemp,yftemp,-z4temp,-z4temp,Axx,Ayy,Azz,intsteps)); %bit that does correct field
  B(counter,(4:6))=B(counter,(4:6))+(I4*linemag(k,xf4,xf4,yftemp,yftemp,-z4temp,-(z4temp+zoff),Axx,Ayy,Azz,intsteps)); % bit that heads out
  B(counter,(4:6))=B(counter,(4:6))+(I4*linemag(k,xs4,xf4,ystemp,ystemp,-(z4temp+zoff),-z4temp,Axx,Ayy,Azz,intsteps)); % bit that heads in
  B(counter,(4:6))=B(counter,(4:6))+(I4*linemag(k,xs4,xf4,yftemp,ystemp,-(z4temp+zoff),-(z4temp+zoff),Axx,Ayy,Azz,intsteps)); % bit does incorrect field
      
   wirecount4=wirecount4+1;
        if (wirecount4==N4)
        break
        end

end  
end  

B(counter,(4:6))=B(counter,(4:6))+[addedbias 0 0];
end
if axiss(ll)==1
Bmag=plotline(B,steps,ll);
hold on
end
counter=0;
end

grid on

