function Bout=Bfield_path_numerical(path_sym,tlim,dt,curr,rot_vec,xyz)
% Bout = Bfield_path_numerical
% B field calculator for arb path of current 
% example
% syms t
% path_sym=symfun([0,0,t],t);
% Bfield_path_numerical(path_sym,tlim,dt,curr,rot_vec,xyz)
%
% Bout is a 3x1 cell-array: {Bx,By,Bz} defined at points x,y,z
% TO DO
%try and speed up compared to the vectorized version
%validate
%warning('this calculation is not validated, do not trust')


% physical constants
%global const
%mu_0=const.mu0;  
%hack to get to work in parfor
mu_0=1.2566370614*10^-6;

dpath_sym=diff(path_sym,'t');
path_fun=matlabFunction(path_sym,'Vars',{'t'});
dpath_fun=matlabFunction(dpath_sym,'Vars',{'t'});

rot_mat=rotationVectorToMatrix(rot_vec);
rev_rot_mat=rotationVectorToMatrix(-rot_vec); %reverse rotation vector
xyz=xyz*rot_mat;

Bout_trapz=xyz*nan; %initalize
tvec=linspace(tlim(1),tlim(2),ceil(range(tlim)/dt))';
%define the parametric curve (with t is the parameter) for the curren
wire_pos=vec_fun(path_fun,tvec); %slow point!
%define the tangent to this curve
dLdt=vec_fun(dpath_fun,tvec); %slow point !

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
