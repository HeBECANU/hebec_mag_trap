clear all
close all

[Hemass,Hegamma,Helambda,Helife,HeIs,Hemu,hbar,kb,Hek,g]=Heconst;

f=560;
ww=f*2*pi;
l=-.01/100:1/100000:.01/100;
BBB=(0.5*Hemass*ww^2*l.^2)/(Hemu*2*pi*hbar)

plot(l*100,BBB,'ok')

db=BBB(1)-BBB(2);
w=sqrt((2*2*pi*Hemu*db*hbar)/(Hemass*abs(l(1)^2-l(2)^2)));

f=w/(2*pi)