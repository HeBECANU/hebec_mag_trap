clear all
%close all

B0=2
Bd=10000;
Bdd=700000;

%x=-.02:.0001:.02;
%y=-.02:.0001:.02;
z=-.12:.0001:.12;
x=0;y=0;
Bx=Bd*x+Bdd/2*(-x.*z);
By=-Bd*y+Bdd/2*(-y.*z);
Bz=B0+Bdd/2*(z.^2-0.5*(x.^2+y.^2));
Bmag=sqrt(Bx.^2+By.^2+Bz.^2);
plot(z,Bmag)