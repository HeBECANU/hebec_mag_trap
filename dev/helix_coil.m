function  B_out=helix_coil(N,r,c,pts);
t = (0:0.1:2*pi*N)';
dl = [-r*sin(t),r*cos(t),c*ones(length(t),1)];
l = [r*cos(t),r*sin(t),c*t];
B_out = biot_savart(l,dl,t,pts);
end
