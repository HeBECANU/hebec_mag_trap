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

%%%%%%%%%%%%%%%%%First coil set %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%z-dist. = 4.5mm (dist from centre to glass face)+ 2.6 (glass thickness)+
%           0.8 (thickness of metal former) +.25( thickness of rim)
%           +0.2(gap between glass and former)
zdd=4.6+2.7+0.8+.2+.25+0e-3;%5+2.6+0.75+0.25+.2;
z1=2*zdd/1000; % distance just before the first wire i.e. ->|o
R1=.75/100;%4.6228
I1=10*(1-0); %20 -> tight, 15 -> high bias for experiment
N1=10;%23.5

%%%%%%%%%%%%%%%%%Second coil set %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
z2=2*zdd/1000;
R2=0.7/100;%4.6228
I2=10*(1+0); %20 -> tight, 30 -> high bias for experiment
N2=19;%23.5
x2=1.85/100; % use 1.95 for higher trap!

%assume the wire coating is around 50 microns thick

%%%%%%%%%%%%%%%%%%%%%%%Work out power consumed


%%%% Note 1T = 10^4 G %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%First CL calc  .%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% over what range? corrdinates in mm
range=5;
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

testpointx=0.46/100;%0.5/100 -> tight 0.23/100 -> high bias experiment
testpointy=0/100;
testpointz=0/100;
addedbias=0e-4;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%What axis do we want plotted%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       [x y z xy]
axiss = [1 1 0 0]

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
        
        if I1~= 0
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%Extended integration for on-axis coil set%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            B(counter,(4:6))=B(counter,(4:6))+(N1*I1*ahmag(k,R1,z1,Axx,Ayy,Azz,intsteps));
        end
        if I2~= 0
            B(counter,(4:6))=B(counter,(4:6))+(N2*I2*ahmagoff(k,x2,R2,z2,Axx,Ayy,Azz,intsteps));
        end
        
        B(counter,(4:6))=B(counter,(4:6))+[-.82e-4 1e-4 0e-4];
    end
    if axiss(ll)==1
        Bmag=plotline(B,steps,ll);
        hold on
        counter=0;
    end
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
    
    f=700;
    w=f*2*pi;
    fp(1)=w;
    fitt=fminsearch('harm_fit',fp,[],ddd,BBB);
    v=(0.5*Hemass*fitt(1)^2*ddd.^2)/(Hemu*2*pi*hbar);
    plot(ddd,BBB,'ok',ddd,v,'r')
    fitt(1)/(2*pi)
end