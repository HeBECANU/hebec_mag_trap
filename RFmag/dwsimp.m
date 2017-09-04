clear all
close all

GtoT=1e-4;
Brf=1.1*GtoT;
G=2000*100*GtoT;
wrf=2*pi*0.5e6;
B1=.75*GtoT;

x =1e-10:.1e-6:20e-6;
y =1e-10:.1e-6:20e-6;
rho=sqrt(x.^2+y.^2);
phi=atan(y./x);
Bsx=G*x*GtoT;
Bsy=G*y*GtoT;
Bm =sqrt(Bsx.^2+Bsy.^2+B1^2);
gf=.5;
mf=2;
mub=9.274e-24;
hbar=1.054e-34;
kb=1.38e-23;
kap =gf*mub;

delta=(kap*B1/hbar)-wrf;
Bc=2*sqrt(B1*hbar*delta/abs(kap))
rho0=1/(sqrt(2)*G)*sqrt(Brf^2-Bc^2)
d=delta/(2*pi)
Vdw=mf*kap*sqrt((Bm-(hbar*wrf/abs(kap))).^2+((Brf./(2*Bm)).^2.*(B1^2+(G^2*rho.^2.*sin(phi).^2))))/kb;
B0=Brf/2*sqrt(1+hbar*delta/(abs(kap)*B1));
Vdw2=mf*kap*sqrt(G^4/B1^2*(((rho.^2-rho0^2)/2).^2)+B0^2)/kb;
plot(x,Vdw,'r',x,Vdw2,'b')