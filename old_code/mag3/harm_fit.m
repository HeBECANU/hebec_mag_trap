function err= harm_fit(p,s,c)

w=p(1);

[Hemass,Hegamma,Helambda,Helife,HeIs,Hemu,hbar,kb,Hek,g]=Heconst;

v=(0.5*Hemass*w^2*s.^2)/(Hemu*2*pi*hbar);
v=v-c;
err=sum(v.^2);