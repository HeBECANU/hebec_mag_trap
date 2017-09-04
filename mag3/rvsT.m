clear all
close all
T=0.5e-6:1e-7:500e-6;

r=0.068*sqrt(T);

plot(T/1e-6,r/1e-3)
ylabel('radius (mm)');
xlabel('Temperature (microK)')

w1f=2*pi*51.9;
w2f=2*pi*74.9;
wbarf=(w1f*w2f^2)^(1/3);
N0=10000;
a=20e-9;
[Hemass,Hegamma,Helambda,Helife,HeIs,Hemu,hbar,kb,Hek,g]=Heconst;
s1f=sqrt(hbar/(Hemass*w1f));
sbarf=sqrt(hbar/(Hemass*wbarf));
x1f=s1f*sqrt(wbarf/w1f)*(15*N0*a/sbarf)^(1/5)

w1i=2*pi*66.7;
w2i=2*pi*1656;
wbari=(w1i*w2i^2)^(1/3);
N0=10000;
a=20e-9;
s1i=sqrt(hbar/(Hemass*w1i));
s2i=sqrt(hbar/(Hemass*w2i));
sbari=sqrt(hbar/(Hemass*wbari));

x1i=s1i*sqrt(wbari/w1i)*(15*N0*a/sbari)^(1/5)
x2i=s2i*sqrt(wbari/w2i)*(15*N0*a/sbari)^(1/5)