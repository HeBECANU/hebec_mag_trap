clear all
B0=3.2;
a=1667;
r=-.05:.0001:.05;
B=sqrt((a*r).^2+B0^2);

plot(r,B,'g*');