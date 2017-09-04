%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%onetrajmag.M%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This program simulates the trajectory of one atom out of the mag guide
% into the magnetic trap, defined by asq (antihelmquic).  JUST uses the
% magnetic forces on an atom assumed to be in the +1 state.
% the long axis of the magnetic trap lies on the x-axis

clear all
close all
tic
%%%%%%%%%%%%%%%%%%%%%%Helium Parameters%%%%%%%%%%%%%%%%%%%%%%%%
[Hemass,Hegamma,Helambda,Helife,HeIs,Hemu,hbar,kb,Hek,g]=Heconst;

%%%%%%%%%%%%%%%%%Simulation Parameters%%%%%%%%%%%%%%%%%%%%%%%%%
T=10e-3; % maximum simulation time
dt =1e-5; % conservative time interval

%%%%%%%%%%%%%%%%%%%Initial conditions%%%%%%%%%%%%%%%%%%%%%%%%%%
vx=25;vy=-1;vz=.9; x=-5/100;y=.2/100;z=.2/100;
options=odeset('RelTol',1e-2,'AbsTol',1e-3,'maxstep',dt);

t=0;
counter=0;
while x<.5/100 & t <T
      %%%initial values for each dt%%%%
    x0=[vx x];
    y0=[vy y]; 
    z0=[vz z];
 
    Tspan=[t,t+dt];
    counter=counter+1;
    
    %%%%%%%%%%%%%%%%%Work out magnetic acceleration%%%%%%%%%%%%%%%%%%%%%%%%%%%
    a=Hemu*hbar*2*pi*gradB(x0(2),y0(2),z0(2),1e-6)/Hemass;
    
          [ti,xx] = ODE45(@newton,Tspan,x0,options,a(1));
          [ti,yy] = ODE45(@newton,Tspan,y0,options,a(2));
          [ti,zz] = ODE45(@newton,Tspan,z0,options,a(3));
        
          t=ti(length(ti));
          vx=xx(length(xx),1);
          x=xx(length(xx),2); 
          vy=yy(length(yy),1);
          y=yy(length(yy),2); 
          vz=zz(length(zz),1);
          z=zz(length(zz),2); 
                  
              velx(counter)=vx;
              posx(counter)=x;
               vely(counter)=vy;
              posy(counter)=y;
               velz(counter)=vz;
              posz(counter)=z;
              time(counter)=t;
end
toc

subplot(3,2,1)
plot(time,velx)
xlabel('time (s)')
Ylabel('x-velocity (m/s) ')

subplot(3,2,2)
plot(time,posx)
xlabel('time (s)')
Ylabel('x-position (m) ')    

subplot(3,2,3)
plot(time,vely)
xlabel('time (s)')
Ylabel('y-velocity (m/s) ')

subplot(3,2,4)
plot(time,posy)
xlabel('time (s)')
Ylabel('y-position (m) ')

subplot(3,2,5)
plot(time,velz)
xlabel('time (s)')
Ylabel('z-velocity (m/s) ')

subplot(3,2,6)
plot(time,posz)
xlabel('time (s)')
Ylabel('z-position (m) ')
    
