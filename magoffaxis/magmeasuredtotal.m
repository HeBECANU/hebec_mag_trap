

x= [0 1 2 3 4 5 6 7 8 9 10 11 12 13 13.5 14 14.5 15 15.5 16 16.5 17 17.5 18 18.5 19 19.5 20 20.5 21 21.5 22 22.5 23 23.5];
BB=[48.34 48.45 47.23 44.69 41 36.39 31.16 25.60 20.02 14.64 9.78 5.68 2.55 0.58 0.09 0.03 0.15 0.68 1.48 2.48 3.70 5.06...
        6.46 7.82 9.15 10.20 11.04 11.58 11.77 11.58 10.99 9.97 8.58 6.74 4.60];

x=(x-14)/10;


hold on
plot (x,BB,'+b')
grid on
