% program to test the relative strength of a cloverleaf trap to a BiQUIC
% assume each coil has only 1 amp-turn

close all
clear all

Mu0 = 4*pi*1e-7;
k=(Mu0)/(4*pi);
intsteps=500;

zdd=4.6+2.7+0.8+.2+.25;%5+2.6+0.75+0.25+.2;
z1=2*zdd/1000; % distance just before the first wire i.e. ->|o

%%%BIQUIC parameters%%%%
R1=.9/100;%4.6228
I1=1;
I2=-1.2;
z2=z1;
x2=.9/100;
R2=2.4/100;
xgi=-1.5/100;xgf=1.5/100;
yg=.5/100;
Ig=0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%lets do BiQUIC first%%%%%%

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
steps=80;
[xvec,yvec,zvec]=setupcoords(xmin,xmax,ymin,ymax,zmin,zmax,steps);
B=zeros(21,6);
counter=0;
inclend=0;


testpointx=-.5/100;%0.5/100 -> tight 0.23/100 -> high bias experiment
testpointy=0/100;
testpointz=0/100;
addedbias=0e-4;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%What axis do we want plotted%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       [x y z xy]
axiss = [0 0 0 0]

for ll=1:length(axiss)
for ttt=1:steps+1
% check if this axis is wanted %
    if axiss(ll)==0
        break
    end
 
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

B(counter,(4:6))=(I1*ahmag(k,R1,z1,Axx,Ayy,Azz,intsteps))+(I2*ahmagoff(k,x2,R2,z2,Axx,Ayy,Azz,intsteps));
B(counter,(4:6))=B(counter,(4:6))+(Ig*linemag(k,xgi,xgf,-yg,-yg,z2/2,z2/2,Axx,Ayy,Azz,intsteps)); 
B(counter,(4:6))=B(counter,(4:6))+(Ig*linemag(k,xgf,xgi,+yg,+yg,z2/2,z2/2,Axx,Ayy,Azz,intsteps)); 
B(counter,(4:6))=B(counter,(4:6))+(Ig*linemag(k,xgf,xgi,-yg,-yg,-z2/2,-z2/2,Axx,Ayy,Azz,intsteps)); 
B(counter,(4:6))=B(counter,(4:6))+(Ig*linemag(k,xgi,xgf,+yg,+yg,-z2/2,-z2/2,Axx,Ayy,Azz,intsteps)); 
if inclend==1
    B(counter,(4:6))=B(counter,(4:6))+(Ig*linemag(k,xgf,xgf,-yg,+yg,z2/2,z2/2,Axx,Ayy,Azz,intsteps)); 
    B(counter,(4:6))=B(counter,(4:6))+(Ig*linemag(k,xgi,xgi,+yg,-yg,z2/2,z2/2,Axx,Ayy,Azz,intsteps)); 
    B(counter,(4:6))=B(counter,(4:6))+(Ig*linemag(k,xgi,xgi,-yg,+yg,-z2/2,-z2/2,Axx,Ayy,Azz,intsteps)); 
    B(counter,(4:6))=B(counter,(4:6))+(Ig*linemag(k,xgf,xgf,+yg,-yg,-z2/2,-z2/2,Axx,Ayy,Azz,intsteps)); 
end
end
if axiss(ll)==1
Bmag=plotline(B,steps,ll);
hold on
end
counter=0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%First CL calc  .%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

zdd=4.6+2.7+0.8+.2+.25;%5+2.6+0.75+0.25+.2;
z1=2*zdd/1000; % distance just before the first wire i.e. ->|o

%%%clover parameters%%%%
Rc=.75/100;%4.6228
Ic=1;
Ib=1;
Rb=2.3/100;
Rcl=.3/100;
Icl=4;
RR=1.25/100;

ox=abs(RR*cos(45/180*pi));
oy=abs(RR*sin(45/180*pi));

% over what range? corrdinates in mm
range=10;
xmax=1*range;
xmin=-1*range;
ymax=1*range;
ymin=-1*range;
zmax=1*range;
zmin=-1*range;
steps=80;
[xvec,yvec,zvec]=setupcoords(xmin,xmax,ymin,ymax,zmin,zmax,steps);
B=zeros(21,6);
counter=0;

testpointx=0/100;%0.5/100 -> tight 0.23/100 -> high bias experiment
testpointy=0/100;
testpointz=0/100;
addedbias=0e-4;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%What axis do we want plotted%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       [x y z xy]
axiss = [1 1 1 0]

for ll=1:length(axiss)
for ttt=1:steps+1
% check if this axis is wanted %
    if axiss(ll)==0
        break
    end
 
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

B(counter,(4:6))=Ic*curvmag(k,Rc,z1,Axx,Ayy,Azz,intsteps);
B(counter,(4:6))=B(counter,(4:6))+(Ib*biasmag(k,Rb,z1,Axx,Ayy,Azz,intsteps));
B(counter,(4:6))=B(counter,(4:6))+(Icl*fourcoiloffsetxandy(k,Rcl,ox,oy,z1,Axx,Ayy,Azz,intsteps));
B(counter,(4:6))=B(counter,(4:6))-(Icl*fourcoiloffsetxandy(k,Rcl,ox,oy,-z1,Axx,Ayy,Azz,intsteps));
end
if axiss(ll)==1
Bmag=plotline(B,steps,ll);
hold on
end
counter=0;
end