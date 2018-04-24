clear all
close all

D=.855/100;

R=0:.0001:4/100;

y=R.^2.*(4*D^2-R.^2)./(D^2+R.^2).^(7/2);

[i j]=max(y)
rr=R(j)

ee=3*R.^4+8*D^4-24*D^2*R.^2;

plot(R,y/max(y),R,ee/max(ee))

[Hemass,Hegamma,Helambda,Helife,HeIs,Hemu,hbar,kb,Hek,g]=Heconst;
mu0 = 4*pi*1e-7;

rr=0.007;
gg=(mu0*1.5*rr^2*(4*D^2-rr^2)/(D^2+rr^2)^(7/2));
gg=gg*250;
w=sqrt(gg*hbar*2*pi*Hemu*1e4*2/Hemass)

f=w/(2*pi)