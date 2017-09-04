clear all
close all

x=-2:.001:2;
Bdash1=100;
B1=Bdash1*x;

B1mag=sqrt(B1.^2)

x2=x-.5;
Bdash2=200
B2=Bdash2*x2;
B2mag=sqrt(B2.^2);

B=B1+B2;
Bmag=sqrt(B.^2);
plot(x,B2mag,'r',x,B1mag,'b',x,Bmag,'g')
axis([-1 1 0 100]);