

x= [-2.5 -2 -1.5 -1 -.5 0 0.5 1 1.5 2 2.5 3 3.5 4 4.5 5 5.5 6 6.5 7 7.5 8 8.5 9 9.5 10 10.5 11 11.5 12 12.5 13 13.5 14 14.5 15 15.5 16 16.5 17 17.5 18 18.5 19 19.5 20 20.5 21 21.5];
BB=[11.3 6.39 1.55 3.22 8.01 12.91 18.04 23.38 28.98 34.83 41.12 47.79 55 62.68 70.73 79.36 88.27 97.62 107.12 116.81 126.36 135.73 144.63 152.93 160.59 167.34 173.29 178.19 181.98...
        184.68 186.3 186.79 186.1 184.33 181.43 177.5 172.5 166.6 159.75 152.28 144.06 135.34 126.48 117.45 108.46 99.66 91.26 83.29 75.8];                                                                                                                                                                                                                                        

x=(x+1.3)/10;

hold on
plot (x,2*BB,'+b')
grid on


x= [0.5 1 1.5 2 2.5 3 3.5 4 4.5 5 5.5 6 6.5 7 7.5 8 8.5 9 9.5 10 10.5 11 11.5 12 12.5 13 13.5 14 14.5 15 15.5 16 16.5 17 17.5 18 18.5 19 19.5 20 20.5 21 21.5];
BB=[39.22 40.77 42.14 43.44 44.61 45.53 46.19 46.52 46.52 46.14 45.35 44.16 42.57 40.5 38.06 35.23 32.0 28.43 24.63 20.47 16.18 6.97 2.26...
        2.53 7.3 12.06 16.75 21.26 25.64 29.74 33.64 37.2 40.46 43.2 45.8 47.86 49.5 50.7 51.4 51.8 51.77 51.4 50.8];

x=(x-11.8)/10;

plot (-1*x,BB,'+b')
grid on
hold off


