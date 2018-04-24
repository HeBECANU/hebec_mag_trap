clear all
close all


Bd=100*100;
Bdd=40*100*100/2;
Bias=1e-4;
ll=3;
if ll==1
x=-1/100:1/100000:1/100;
y=0;
z=0;
end
if ll==2
y=-1/100:1/100000:1/100;
x=0;
z=0;
end
if ll==3
z=-1/100:1/100000:1/100;
y=0;
x=0;
end

Bx0=0;
By0=5;
Bz0=.5;
    
Bx=Bx0+Bd*x-Bdd*x.*z;
By=By0-Bd*y-Bdd*y.*z;
Bz=Bz0+Bdd*(z.^2-(x.^2+y.^2)/2)

B=sqrt(Bx.^2+By.^2+Bz.^2);

if ll==1
plot(x,B)
end
if ll==2
plot(y,B)
end

if ll==3
plot(z,B)
end

