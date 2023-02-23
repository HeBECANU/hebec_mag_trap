function Bout=Bfield_path_numeric_order_1(path_sym,tlim,dt,curr,rot_vec,xyz)
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
% for parfor usage we specify the physical constants
mu_0=1.2566370614*10^-6;

dpath_sym=diff(path_sym,'t');
path_fun=matlabFunction(path_sym,'Vars',{'t'});

rot_mat=rotationVectorToMatrix(rot_vec);
rev_rot_mat=rotationVectorToMatrix(-rot_vec); %reverse rotation vector
xyz=xyz*rot_mat;

Bout_sum=xyz*nan; %initalize
tvec=linspace(tlim(1),tlim(2),ceil(range(tlim)/dt))';
%define the parametric curve (with t is the parameter) for the curren
wire_pos=vec_fun(path_fun,tvec); %slow point!

for ii=1:size(xyz,1)
     %difference between interogation point current stick start
    bvec=xyz(ii,:)-wire_pos(1:end-1,:);
    %difference between interogation point current stick end
    cvec=xyz(ii,:)-wire_pos(2:end,:);
    %vector of the differnece between start and end
    avec=wire_pos(2:end,:)-wire_pos(1:end-1,:);
    %calculate the magnetic feild from each current stick
    % from eq 22 we caulate the feild from each stick
    c_cross_a=cross(cvec,avec,2);
    Bvec=(c_cross_a./(vecnorm(c_cross_a,2,2).^2)) .*( ( dot(avec,cvec,2)./vecnorm(cvec,2,2) ) - ( dot(avec,bvec,2)./vecnorm(bvec,2,2) ) ) ;
    %combine them all
    Bout_sum(ii,:)=sum(Bvec,1);
end
Bout=((mu_0*curr)/(4*pi)).*Bout_sum;

Bout=Bout*rev_rot_mat;


end
