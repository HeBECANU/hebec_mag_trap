clear all
close all

%%%results for Ioffe coil

II = [1 2 3 4 5 10 12 14 16 18 20 22 23 24 25 26 27 28 29 30 31 33 35 37 40];
VI= [70e-3 140e-3 208e-3 280e-3 350e-3 .703 .851 .999 1.161 1.320 1.484 1.66 1.747 1.84 1.934 2.030 2.127 ...
        2.224 2.328 2.440 2.551 2.768 2.994 3.233 3.6];

PI= VI.*II;

RI=VI./II;

RrtI=.07;
TI=(RI/RrtI-1)/.004;

IIt=[1 10 16 20 25 30 40];
VIt=[70e-3 .703 1.161 1.484 1.934 2.440 3.6];
Ti= [19 22.5 28.5 33.8 43 53 77]-19;
PIt=IIt.*VIt;

plot(PI,TI,'.b')
hold on
plot(PIt,Ti,'*r')
hold off
ylabel('Temp increase (degrees)')
xlabel('Power (Watts)')

pause

IQ= [10 12 14 16 20 22 25 30 32 35 38 40];
VQ=[.438 .529 .621 .713 .906 1.010 1.168 1.448 1.570 1.764 1.968 2.121];
PQ=IQ.*VQ;

RQ=VQ./IQ;

RrtQ=0.043;
TQ=(RQ/RrtQ-1)/.004;

Tq=[20.7 21.3 22.1 23 25.1 26.2 28.5 32.3 34.3 37.5 41.1 43.7]-19;
plot(PQ,TQ,'.b')
hold on
plot(PQ,Tq,'*r')
hold off

ylabel('Temp increase (degrees)')
xlabel('Power (Watts)')
