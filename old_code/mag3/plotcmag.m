clear all


load xout.dat
load yout.dat
load zout.dat
load xyout.dat



plot(xout(:,1),xout(:,5),'r',yout(:,1),yout(:,5),'b',xyout(:,1),xyout(:,5),'.g',zout(:,1),zout(:,5),'y');
grid