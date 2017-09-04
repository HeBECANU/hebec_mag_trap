%a program to calculate the field for a 
%QUIC trap. The trap comprises two sets of anti-helmholtz coils
% offset from each other
% This is an extended coil model

%same as antihelmquic, but now we have added the lead in and lead out wires
%  We assume the lead in and out wires are collinear

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
intsteps=500;

alpha=17e-6; % linear thermal expansion coefficient/ degree, for copper
dT=0;%ten degree change!!
dx=alpha*dT*(0.5/1000);
%%%%%%%%%%%%%%%%%First coil set %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
z1=2*.85/100; % distance to the just before the first wire i.e. ->|o
R1=.75/100;%4.6228
I1=20; %20 -> tight, 15 -> high bias for experiment
N1=10;%23.5
thick1=.5/1000+dx; % wire thickness 
numz1s=10;
numr1s=1;
%%%%%%%%%%%%%%%%First coil lead wires%%%%%%%%%%%%%%%%%%%%%%
thickw=.5/1000+dx; % wire thickness 
w1startx =0; w1starty=R1+thick1+thickw/2;; w1startz=z1;
w1stopx =1e-12; w1stopy=R1+thick1+thickw/2+1e-12;; w1stopz=z1+(numz1s*thick1);

%%%%%%%%%%%%%%%%%Second coil set %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
z2=2*0.85/100;
R2=0.7/100;%4.6228
I2=20; %20 -> tight, 30 -> high bias for experiment
N2=17  ;%23.5
x2=1.85/100;
thick2=.5/1000+dx;
numz2s=17;
numr2s=1; 
%%%%%%%%%%%%%%%%Second coil lead wires%%%%%%%%%%%%%%%%%%%%%%
thickw=.5/1000+dx; % wire thickness 
w2startx =0; w2starty=R2+thick2+thickw/2;; w2startz=z2;
w2stopx =1e-12; w2stopy=R2+thick2+thickw/2+1e-12; w2stopz=z2+(numz2s*thick2);
%assume the wire coating is around 50 microns thick

%%%%%%%%%%%%%%%%%%%%%%%Work out power consumed
rho=1.75e-8; %ohms m
wthick1=(thick1/2)-50e-6;
a1=pi*wthick1^2;
l1=2*N1*2*pi*R1;

wthick2=(thick2/2)-50e-6;
a2=pi*wthick2^2;
l2=2*N2*2*pi*R2;

Res1=rho*l1/a1
Res2=rho*l2/a2

Res=Res1+Res2
V=I1*Res


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
steps=80;
[xvec,yvec,zvec]=setupcoords(xmin,xmax,ymin,ymax,zmin,zmax,steps);
counter=0;
wirecount=0;
B=zeros(21,6);

testpoint=0.49/100;%0.5/100 -> tight 0.23/100 -> high bias experiment
onlyx=0;
addedbias=0e-4;
if onlyx ==0
% first calculate for z-axis
for ttt=1:steps+1

    wirecount1=0;
    wirecount2=0;
% assign test points
    Axx=testpoint;
    Ayy=0;
    Azz=zvec(ttt);
    counter=counter+1;

B(counter,1)=Axx;
B(counter,2)=Ayy;
B(counter,3)=Azz;
B(counter,(4:6))=zeros(1,3);

%%%%%Lead in wire for first coil%%%%%%%%%
B(counter,(4:6))=B(counter,(4:6))+(I1*linemag(k,w1startx,w1stopx,w1starty,w1stopy,w1startz,w1stopz,Axx,Ayy,Azz,intsteps));

%%%%Extended integration for on-axis coil%%%%%%%%%%%%%%%%%%%%%%%%
for ss=1:numz1s
    if (wirecount1==N1)
    break
    end
    z1temp=z1+(thick1/2)+(thick1*(ss-1)); % distance to the centre of each wire 
for tt =1:numr1s
    R1temp=R1+(thick1/2)+(thick1*(tt-1));

    B(counter,(4:6))=B(counter,(4:6))+(I1*ahmag(k,R1temp,z1temp,Axx,Ayy,Azz,intsteps));
    wirecount1=wirecount1+1;
        if (wirecount1==N1)
        break
        end
end
end
%%%%%Lead in wire for second coil%%%%%%%%%
B(counter,(4:6))=B(counter,(4:6))+(I2*linemag(k,w2startx,w2stopx,w2starty,w2stopy,w2startz,w2stopz,Axx,Ayy,Azz,intsteps));
%%%%Extended integration for off-axis coil%%%%%%%%%%%%%%%%%%%%%%%%
for ss=1:numz2s
    if (wirecount2==N2)
    break
    end
    z2temp=z2+(thick2/2)+(thick2*(ss-1));   
    
    for tt =1:numr2s
    R2temp=R2+(thick2/2)+(thick2*(tt-1));  

        B(counter,(4:6))=B(counter,(4:6))+(I2*ahmagoff(k,x2,R2temp,z2temp,Axx,Ayy,Azz,intsteps));
        wirecount2=wirecount2+1;
        if (wirecount2==N2)
        break
        end
    end
end  
B(counter,(4:6))=B(counter,(4:6))+[addedbias 0 0];
end
Bmag=plotline(B,steps,3);
hold on
counter=0;

% now calculate for y-axis
for ttt=1:steps+1

    wirecount1=0;
    wirecount2=0;

    % assign test points
    Axx=testpoint;
    Ayy=yvec(ttt);
    Azz=0;

    counter=counter+1;
    B(counter,1)=Axx;
    B(counter,2)=Ayy;
    B(counter,3)=Azz;
    B(counter,(4:6))=zeros(1,3);
%%%%%Lead in wire for first coil%%%%%%%%%
B(counter,(4:6))=B(counter,(4:6))+(I1*linemag(k,w1startx,w1stopx,w1starty,w1stopy,w1startz,w1stopz,Axx,Ayy,Azz,intsteps));
for ss=1:numz1s
    if (wirecount1==N1)
    break
    end
    z1temp=z1+(thick1/2)+(thick1*(ss-1));  
for tt =1:numr1s
    R1temp=R1+(thick1/2)+(thick1*(tt-1));

    B(counter,(4:6))=B(counter,(4:6))+(I1*ahmag(k,R1temp,z1temp,Axx,Ayy,Azz,intsteps));
    wirecount1=wirecount1+1;
    if (wirecount1==N1)
    break
    end
end
end
%%%%%Lead in wire for second coil%%%%%%%%%
B(counter,(4:6))=B(counter,(4:6))+(I2*linemag(k,w2startx,w2stopx,w2starty,w2stopy,w2startz,w2stopz,Axx,Ayy,Azz,intsteps));
for ss=1:numz2s
    if (wirecount2==N2)
    break
    end
    z2temp=z2+(thick2/2)+(thick2*(ss-1));   
    for tt =1:numr2s
    R2temp=R2+(thick2/2)+(thick2*(tt-1));  

    B(counter,(4:6))=B(counter,(4:6))+(I2*ahmagoff(k,x2,R2temp,z2temp,Axx,Ayy,Azz,intsteps));
    wirecount2=wirecount2+1;
    if (wirecount2==N2)
    break
    end
    end
end
B(counter,(4:6))=B(counter,(4:6))+[addedbias 0 0];
end

Bmag=plotline(B,steps,2);
counter=0;

end

% now calculate for x-axis
for ttt=1:steps+1

wirecount1=0;
wirecount2=0;
% assign test points
    Axx=xvec(ttt)+testpoint;
    Ayy=0;
    Azz=0;
    counter=counter+1;

B(counter,1)=Axx-testpoint;
B(counter,2)=Ayy;
B(counter,3)=Azz;
B(counter,(4:6))=zeros(1,3);
%%%%%Lead in wire for first coil%%%%%%%%%
B(counter,(4:6))=B(counter,(4:6))+(I1*linemag(k,w1startx,w1stopx,w1starty,w1stopy,w1startz,w1stopz,Axx,Ayy,Azz,intsteps));
for ss=1:numz1s
    if (wirecount1==N1)
    break
    end
    z1temp=z1+(thick1/2)+(thick1*(ss-1));  
for tt =1:numr1s
    R1temp=R1+(thick1/2)+(thick1*(tt-1));

        B(counter,(4:6))=B(counter,(4:6))+(I1*ahmag(k,R1temp,z1temp,Axx,Ayy,Azz,intsteps));
    wirecount1=wirecount1+1;
    if (wirecount1==N1)
    break
    end
end
end
%%%%%Lead in wire for second coil%%%%%%%%%
B(counter,(4:6))=B(counter,(4:6))+(I2*linemag(k,w2startx,w2stopx,w2starty,w2stopy,w2startz,w2stopz,Axx,Ayy,Azz,intsteps));
for ss=1:numz2s
    if (wirecount2==N2)
    break
    end
    z2temp=z2+(thick2/2)+(thick2*(ss-1));   
for tt =1:numr2s
    R2temp=R2+(thick2/2)+(thick2*(tt-1));  

    B(counter,(4:6))=B(counter,(4:6))+(I2*ahmagoff(k,x2,R2temp,z2temp,Axx,Ayy,Azz,intsteps));
    wirecount2=wirecount2+1;
    if (wirecount2==N2)
    break
    end
end
end
B(counter,(4:6))=B(counter,(4:6))+[addedbias 0 0];
end

Bmag=plotline(B,steps,1);
counter=0;


if onlyx==0
% now calculate for xy-axis
for ttt=1:steps+1   
    wirecount1=0;
    wirecount2=0;

% assign test points
    Axx=xvec(ttt)+testpoint;
    Ayy=yvec(ttt);
    Azz=0;
    counter=counter+1;
    B(counter,1)=Axx-testpoint;
    B(counter,2)=Ayy;
    B(counter,3)=Azz;
    B(counter,(4:6))=zeros(1,3);
    %%%%%Lead in wire for first coil%%%%%%%%%
B(counter,(4:6))=B(counter,(4:6))+(I1*linemag(k,w1startx,w1stopx,w1starty,w1stopy,w1startz,w1stopz,Axx,Ayy,Azz,intsteps));
for ss=1:numz1s
    if (wirecount1==N1)
    break
    end
    z1temp=z1+(thick1/2)+(thick1*(ss-1));  
    for tt =1:numr1s
    R1temp=R1+(thick1/2)+(thick1*(tt-1));

    B(counter,(4:6))=B(counter,(4:6))+(I1*ahmag(k,R1temp,z1temp,Axx,Ayy,Azz,intsteps));
    wirecount1=wirecount1+1;
    if (wirecount1==N1)
    break
    end
    end
end
%%%%%Lead in wire for second coil%%%%%%%%%
B(counter,(4:6))=B(counter,(4:6))+(I2*linemag(k,w2startx,w2stopx,w2starty,w2stopy,w2startz,w2stopz,Axx,Ayy,Azz,intsteps));
for ss=1:numz2s
    if (wirecount2==N2)
    break
    end
    z2temp=z2+(thick2/2)+(thick2*(ss-1));   
    for tt =1:numr2s
    R2temp=R2+(thick2/2)+(thick2*(tt-1));  

    B(counter,(4:6))=B(counter,(4:6))+(I2*ahmagoff(k,x2,R2temp,z2temp,Axx,Ayy,Azz,intsteps));
    wirecount2=wirecount2+1;
    if (wirecount2==N2)
    break
    end
end
end
B(counter,(4:6))=B(counter,(4:6))+[addedbias 0 0];
end   
end

Bmag=plotline(B,steps,4);
grid on
hold off
toc