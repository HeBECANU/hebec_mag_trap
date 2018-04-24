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
%           1 (thickness of plastic former) +.25( thickness of rim)
%           +0.2(gap between glass and former)
zdd=4.5+2.6+1+.2+.25;%5+2.6+0.75+0.25+.2;
z1=2*zdd/1000; % distance just before the first wire i.e. ->|o
R1=.75/100;%4.6228
I1=30*(1-0); %20 -> tight, 15 -> high bias for experiment
N1=10;%23.5
thick1z=0.56/1000+0/1000; % wire thickness in the z direction
thick1r=0.56/1000 %% wire thickness in the r direction
stepsr =1; % integration steps in the r-direction
dr1=thick1r/stepsr; %infintessimal in the r dir.
numz1s=10;

%%%%%%%%%%%%%%%%%Second coil set %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
z2=2*zdd/1000;
R2=0.7/100;%4.6228
I2=-30*(1+0); %20 -> tight, 30 -> high bias for experiment
N2=18 ;%23.5
x2=1.85/100; % use 1.95 for higher trap!
thick2z=0.56/1000+0/1000; % wire thickness in the z direction
thick2r=0.56/1000 %% wire thickness in the r direction
stepsr =1; % integration steps in the r-direction
dr2=thick2r/stepsr; %infintessimal in the r dir.
numz2s=18;


%assume the wire coating is around 50 microns thick

%%%%%%%%%%%%%%%%%%%%%%%Work out power consumed


%%%% Note 1T = 10^4 G %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%First CL calc  .%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% over what range? corrdinates in mm
range=10;
xmax=1*range;
xmin=-1*range;
ymax=1*range;
ymin=-1*range;
zmax=1*range;
zmin=-1*range;
steps=100;
[xvec,yvec,zvec]=setupcoords(xmin,xmax,ymin,ymax,zmin,zmax,steps);
counter=0;
wirecount=0;
B=zeros(21,6);
addwires=0;

testpointx=(0.5)/100;%0.5/100 -> tight 0.23/100 -> high bias experiment
testpointy=0/100;
testpointz=0;%.005/100-3e-4/100;
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
% assign test points
if ll == 1
    Axx=xvec(ttt)+testpointx; Ayy=testpointy; Azz=testpointz;
end
if ll == 2
    Axx=testpointx; Ayy=yvec(ttt)+testpointy; Azz=testpointz;
end
if ll == 3
    Axx=testpointx; Ayy=testpointy; Azz=zvec(ttt)+testpointz;
end
if ll == 4
    Axx=xvec(ttt)+testpointx; Ayy=yvec(ttt)+testpointy; Azz=testpointz;
end
    counter=counter+1;

   B(counter,1)=Axx-testpointx;
    B(counter,2)=Ayy-testpointy;
    B(counter,3)=Azz-testpointz;
    B(counter,(4:6))=zeros(1,3); 

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
B(counter,(4:6))=B(counter,(4:6))+[0e-4 0e-4 0e-4];
end
if axiss(ll)==1
Bmag=plotline(B,steps,ll);
hold on
end
counter=0;
end
[Hemass,Hegamma,Helambda,Helife,HeIs,Hemu,hbar,kb,Hek,g]=Heconst;

plharm=0;
if plharm==1 % plots a comparison harmonic profile
l=-.02/100:1/100000:.02/100;
l=l-(.004/100);
f=70;
ww=f*2*pi;
BBB=(0.5*Hemass*ww^2*l.^2)/(Hemu*2*pi*hbar)+1.26;
hold on
l=l+(.005/100);
plot(l*100,BBB,'ok')
hold off
end
grid on
fittt=0;
if fittt==1 % tries to fit for the harmonic frequency
    fitrange=6;
  
    [i,j]=min(Bmag); % find  minimum
    BBB=Bmag-i;
    if axiss(1)==1
        ddd=B(j-fitrange:j+fitrange,1)-B(j,1);
    end
    if axiss(2)==1
        ddd=B(j-fitrange:j+fitrange,2)-B(j,2);
    end
    if axiss(3)==1
        ddd=B(j-fitrange:j+fitrange,3)-B(j,3);
    end
    
    close all
    BBB=BBB(j-fitrange:j+fitrange);
 
    plot(ddd,BBB,'ok')
    pause
    
    f=200;
   w=f*2*pi;
fp(1)=w;
fitt=fminsearch('harm_fit',fp,[],ddd,BBB); 
    v=(0.5*Hemass*fitt(1)^2*ddd.^2)/(Hemu*2*pi*hbar);
     plot(ddd,BBB,'ok',ddd,v,'r')
     fitt(1)/(2*pi)
end

    
    
    
    
    