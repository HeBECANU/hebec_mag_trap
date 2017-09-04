%%%%%%% subroutine for calculating magnetic field components at (Ax,Ay,Az)
% for all four bars of the magnetic guide.

function[Btotal]=Loadmagend(k,x1start,y1start,z1start,x1stop,y1stop,z1stop,...
        x2start,y2start,z2start,x2stop,y2stop,z2stop,Ax,Ay,Az,intsteps)

% we assumme the user has input the start and stop values so that
% the direction of the current flows from the start to the stop coordinate.


Bline1=maglinemat(k,x1start,y1start,z1start,x1stop,y1stop,z1stop,Ax,Ay,Az);
Bline2=maglinemat(k,x2start,y2start,z2start,x2stop,y2stop,z2stop,Ax,Ay,Az);%   1+    2-


% now lets add all the fields up
Btotal=Bline1-Bline2;

