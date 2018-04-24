clear all
close all

G=100*100; % gradient of the IP field in G/m;
range=50e-6;
counterx=0;
countery=0;
for x=-range:range/20:range
    counterx=counterx+1;
    for y=-range:range/20:range
        countery=countery+1;
Bsx=(G*x);
Bsy=(G*y);
B1=.3;
Ba=.1;
Bb=.1;
Bm=Ba^2+Bb^2;
wrf=2*pi*.5e6;
mub =2*pi*1.4e6; % Rad/G

D1=sqrt(Bsx^2+Bsy^2+B1^2);%in Gauss
gf=0.5;
delta=0;
D2=wrf/(gf*mub); % in Gauss

D=D1-D2;
rho=sqrt(x^2+y^2);
phi=atan(y/x);
alpha=atan(Bb/Ba);
gam=-(gf/abs(gf))*delta;
K=(2*B1*Bm/(8*D1^2));
R1=(B1+(D1*sin(2*alpha)*sin(gam)));
R2=(G^2*rho^2*(1-cos(2*alpha)*cos(2*phi)+sin(2*alpha)*sin(2*phi)*cos(gam)));
Rabsq=K*(R1+R2);
        
mf=2;
V(counterx,countery)=mf*gf*mub*sqrt(D^2+Rabsq);
xx(counterx,countery)=x;
yy(counterx,countery)=y;
end
end


contour(xx,yy,V)
