function  Bout=Bfield_helix(N,r,c,I,pts)
t = (0:0.1:2*pi*N)';
dl = [-r*sin(t),r*cos(t),c*ones(length(t),1)];
l = [r*cos(t),r*sin(t),c*t];
%hack to get to work in parfor
mu_0=1.2566370614*10^-6;
Bout = biot_savart(l,dl,t,pts);
Bout = mu_0*I/(4*pi).*Bout;
end
