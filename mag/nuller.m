clear all
close all

load nulledaxial
plot(aa,bb)

load unnulledaxial
plot(aa,bb,cc,dd)

load nulledupdown
plot(ee,ff)

load unnulledupdown
plot(ee,ff,gg,hh)

load nulledbeam
plot(ii,jj)

load unnulledbeam
plot(ii,jj,kk,ll)

subplot(3,1,1)
plot(aa,bb*500*.009,cc,dd*500*.009)
ylabel('Amplitude (mG)')
xlabel('Time (s)')
subplot(3,1,2)
plot(ee,ff*1000*.009,gg,hh*1000*.009)
ylabel('Amplitude (mG)')
xlabel('Time (s)')
subplot(3,1,3)
plot(ii,jj*1000*.009,kk,ll*1000*.009)
ylabel('Amplitude (mG)')
xlabel('Time (s)')