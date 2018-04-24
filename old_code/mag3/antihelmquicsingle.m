%a program to calculate the field for a 
%QUIC trap. The trap comprises two sets of anti-helmholtz coils
% offset from each other
% This is an extended coil model

%same as antihelmquic, except for a single point (x,y,z)
%returns Bx,By,Bz
function[B]=aqs(x,y,z)


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
z1=2*.9/100; % distance to the just before the first wire i.e. ->|o
R1=.75/100;%4.6228
I1=20; %20 -> tight, 15 -> high bias for experiment
N1=10;%23.5
thick1=.5/1000; % wire thickness 
numz1s=10;
numr1s=1;
%%%%%%%%%%%%%%%%%Second coil set %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
z2=2*0.9/100;
R2=0.7/100;%4.6228
I2=20; %20 -> tight, 30 -> high bias for experiment
N2=16  ;%23.5
x2=1.85/100;
thick2=.5/1000;
numz2s=16;
numr2s=1; 

%%%% Note 1T = 10^4 G %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%First CL calc  .%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    wirecount1=0;
    wirecount2=0;
% assign test points
    Axx=0;
    Ayy=0;
    Azz=0;

B(1)=Axx;
B(2)=Ayy;
B(3)=Azz;
B(4:6)=zeros(1,3);
%%%%Extended integration for on-axis coil%%%%%%%%%%%%%%%%%%%%%%%%
for ss=1:numz1s
    if (wirecount1==N1)
    break
    end
    z1temp=z1+(thick1/2)+(thick1*(ss-1)); % distance to the centre of each wire 
for tt =1:numr1s
    R1temp=R1+(thick1/2)+(thick1*(tt-1));

    B(4:6)=B(4:6)+(I1*ahmag(k,R1temp,z1temp,Axx,Ayy,Azz,intsteps));
    wirecount1=wirecount1+1;
        if (wirecount1==N1)
        break
        end
end
end
%%%%Extended integration for off-axis coil%%%%%%%%%%%%%%%%%%%%%%%%
for ss=1:numz2s
    if (wirecount2==N2)
    break
    end
    z2temp=z2+(thick2/2)+(thick2*(ss-1));   
    
    for tt =1:numr2s
    R2temp=R2+(thick2/2)+(thick2*(tt-1));  

        B(4:6)=B(4:6)+(I2*ahmagoff(k,x2,R2temp,z2temp,Axx,Ayy,Azz,intsteps));
        wirecount2=wirecount2+1;
        if (wirecount2==N2)
        break
        end
    end
end  

B

