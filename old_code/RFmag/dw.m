clear all
close all

G=100*100; % gradient of the IP field in G/m;
range=40e-6;

x=-range:range/1000:range;
y=-range:range/1000:range;
Bsx=(G*x);
Bsy=(G*y);
B1=.5;
Ba=0;
Bb=1.2;
Bm=Ba^2+Bb^2;
wrf=2*pi*.3e6;
mub =2*pi*1.4e6; % Rad/G

D1=sqrt(Bsx.^2+Bsy.^2+B1^2);%in Gauss
gf=0.5;
delta=0;
D2=wrf/(gf*mub); % in Gauss

D=D1-D2;
rho=sqrt(x.^2+y.^2);
phi=atan(y./x);
alpha=atan(Bb/Ba);
gam=-(gf/abs(gf))*delta;
K=(2*B1*Bm./(8*D1.^2));
R1=(B1+(D1*sin(2*alpha)*sin(gam)));
R2=(G^2*rho.^2.*(1-cos(2*alpha)*cos(2*phi)+sin(2*alpha)*sin(2*phi)*cos(gam)));
Rabsq=K.*(R1+R2);
        
mf=2;
V=mf*gf*mub*sqrt(D.^2+Rabsq);

plot(rho,V)

rho=0:50e-6/1000:70e-6;
kap=gf*mub;
del=kap*B1-wrf;
Bc=2*sqrt(B1*del/kap);
rho0=1/(sqrt(2)*G)*sqrt(Bb^2-Bc^2)
B0=Bb/2*sqrt(1+del/(kap*B1));

Vdw=mf*kap*sqrt((G^4/B1^2*((rho.^2-rho0^2)/2).^2)+B0^2);

pause
plot(rho,Vdw)
