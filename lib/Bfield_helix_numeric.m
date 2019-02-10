function Bout=Bfield_helix_numerical(radius,pitch,turn,dt,curr,rot_vec,xyz)
% Bout = Bfield_helix_numerical
% B field calculator for a helix of current 
% Bout is a 3x1 cell-array: {Bx,By,Bz} defined at points x,y,z
% 20190210 validated against Bfield_path_numerical with check_helix_equals_loop.m
% TO DO


% physical constants
%global const
%mu_0=const.mu0;  
%hack to get to work in parfor
mu_0=1.2566370614*10^-6;

rot_mat=rotationVectorToMatrix(rot_vec);
rev_rot_mat=rotationVectorToMatrix(-rot_vec); %reverse rotation vector
xyz=xyz*rot_mat;

Bout_trapz=xyz*nan; %initalize
tlim=2*pi*turn;
tvec=linspace(tlim(1),tlim(2),ceil(range(tlim)/dt))';
%define the parametric curve (with t is the parameter) for the curren
wire_pos=[radius*cos(tvec), radius*sin(tvec), pitch*tvec/(2*pi)]; 
%define the tangent to this curve
dLdt=[-radius*sin(tvec), radius*cos(tvec), tvec*0+pitch/(2*pi)]; 

for ii=1:size(xyz,1)
    %difference between interogation point and the wire position
    Rvec=xyz(ii,:)-wire_pos;
    %calculate the infentesimal b field vector
    dBvec=cross(dLdt,Rvec)./(vecnorm(Rvec,2,2).^3);
    %combine them all
    Bout_trapz(ii,:)=trapz(tvec,dBvec);
end
Bout=((mu_0*curr)/(4*pi)).*Bout_trapz;

Bout=Bout*rev_rot_mat;


end
