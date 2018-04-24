function[f,w]=freqcalc(db,r1,r2)

%%%%%%%%%%%%%%%%%%%%%%Helium Parameters%%%%%%%%%%%%%%%%%%%%%%%%
[Hemass,Hegamma,Helambda,Helife,HeIs,Hemu,hbar,kb,Hek,g]=Heconst;

[Limass,Ligamma,Lilambda,Lilife,LiIs,Limu,hbar,kb,Lik,g]=Liconst;
Hemu
w=sqrt((2*2*pi*Hemu*db*hbar)/(Hemass*abs(r1^2-r2^2)))
f=w/(2*pi)