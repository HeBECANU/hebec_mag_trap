clear all

[Hemass,Hegamma,Helambda,Helife,HeIs,Hemu,hbar,kb,Hek,g]=Heconst;

f=72;
w=2*pi*f;

r=-.5/100:.01/100:.5/100;
B=.5*Hemass*w^2*r.^2/(2*pi*Hemu*hbar);
B0=40;
B=B0+B;
plot(r*100,B)
