%a program to calculate the field for a 
%QUIC trap. The trap comprises two sets of anti-helmholtz coils
% offset from each other
% This is an extended coil model
% this is the programs functional form, so that it can be called
% to calculate the GRAD of the field at a point x,y,z
% used to calculated atomic trajectories
function[dBx,dBy,dBz]=antihelmquicfun(Axx,Ayy,Azz,DD) 

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
z1=2*.85/100; % distance to the just before the first wire i.e. ->|o
R1=.85/100;%4.6228
I1=20; %20 -> tight, 15 -> high bias for experiment
N1=10;%23.5
thick1=0.5/1000; % wire thickness 
numz1s=10;
numr1s=1;
%%%%%%%%%%%%%%%%%Second coil set %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
z2=2*0.85/100;
R2=.6/100;%4.6228
I2=20; %20 -> tight, 30 -> high bias for experiment
N2=22  ;%23.5
x2=1.9/100;
thick2=0.5/1000;
numz2s=11;
numr2s=2; 

%%%% Note 1T = 10^4 G %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%First CL calc  .%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addedbias=0e-4;
%work out gradient in the x-direction
for ttt=1:2
    wirecount1=0;
    wirecount2=0;
    %%%%%%%%STEP either side of Ax to find gradient%%%
if ttt ==1 
    Axxtemp=Axx-DD;
else
    Axxtemp=Axx+DD;
end
B(ttt,1:3)=0;
%%%%Extended integration for on-axis coil%%%%%%%%%%%%%%%%%%%%%%%%
for ss=1:numz1s
    if (wirecount1==N1)
    break
    end
    z1temp=z1+(thick1/2)+(thick1*(ss-1)); % distance to the centre of each wire 
for tt =1:numr1s
    R1temp=R1+(thick1/2)+(thick1*(tt-1));

    B(ttt,1:3)=B(ttt,1:3)+(I1*ahmag(k,R1temp,z1temp,Axxtemp,Ayy,Azz,intsteps));
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

        B(ttt,1:3)=B(ttt,1:3)+(I2*ahmagoff(k,x2,R2temp,z2temp,Axx,Ayy,Azz,intsteps));
        wirecount2=wirecount2+1;
        if (wirecount2==N2)
        break
        end
    end
end  

B(ttt,1:3)=B(ttt,1:3)+[addedbias 0 0];

Bmag(ttt)= sqrt(sum(B(ttt,(1:3)).^2,2));
Bmag(ttt)=Bmag(ttt)/1e-4;
end
Bmag(1)
Bmag(2)
dBx =(Bmag(2)-Bmag(1))/(Axx+DD-(Axx-DD));
dBy=0;dBz=0;





