function Bout=Bfield_loop(radius,curr,rot_vec,xyz)
% Bout = Bfield_coil(R, I, x, y, z)
% B field calculator for single coil (located at origin, pointing Z-axis)
% Bout is a 3x1 cell-array: {Bx,By,Bz} defined at points x,y,z
% TO DO
%CHECK ROTATION OF B VECTOR IS THE RIGHT SIGN
% fix up nan output when on axis of the coil

% physical constants
%global const
%mu_0=const.mu0;  

%hack to get to work in parfor
mu_0=1.2566370614*10^-6;

rot_mat=rotationVectorToMatrix(rot_vec);
rev_rot_mat=rotationVectorToMatrix(-rot_vec); %reverse rotation vector
xyz=xyz*rot_mat;
% cartesian to cylindrical
% theta-rad-z
[theta,rho,z]=cart2pol(xyz(:,1),xyz(:,2),xyz(:,3));

% evaluate the analytic magnetic field solution
% https://doi.org/10.1103/PhysRevA.35.1535
ell_k2=4*radius*rho./((radius+rho).^2+z.^2);      % k^2 elliptic parameter
[K,E]=ellipke(ell_k2);      % complete elliptic integrals of 1st,2nd orders

Bzz=(mu_0*curr./(2*pi)).*(1./sqrt((radius+rho).^2+z.^2)).*(K+E.*(radius^2-rho.^2-z.^2)./((radius-rho).^2+z.^2));
Brr=(mu_0*curr./(2*pi*rho)).*(z./sqrt((radius+rho).^2+z.^2)).*(-K+E.*(radius^2+rho.^2+z.^2)./((radius-rho).^2+z.^2));
Brr(rho==0)=0; %handle the case when the point is on the axis of the loop
Bzz(~isfinite(Bzz))=NaN;    % Inf --> NaN
Brr(~isfinite(Brr))=NaN;

% Reverse transform cyl to original Cart coord (trap centered ref)

[Bout(:,1),Bout(:,2),Bout(:,3)]=pol2cart(theta,Brr,Bzz);
Bout=Bout*rev_rot_mat;


end