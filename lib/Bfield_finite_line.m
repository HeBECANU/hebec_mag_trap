function Bout=Bfield_finite_line(L,I,rot_vec,xyz)
% Bout = Bfield_finite_line(R, I, x, y, z)
% B field calculator for single line of current (located at origin, pointing Z-axis)
% Bout is a 3x1 cell-array: {Bx,By,Bz} defined at points x,y,z
% TO DO
%CHECK ROTATION OF B VECTOR IS THE RIGHT SIGN

%http://web.mit.edu/viz/EM/visualizations/coursenotes/modules/guide09.pdf
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
%see mag_field_of_finite_line mathematica notebook
B_perp=(mu_0*I./(4*pi*rho)).*((L./sqrt(rho.^2+(L-z).^2))...
                             -(z./sqrt(rho.^2+(L-z).^2))...
                             +(z./sqrt(rho.^2+z.^2)));
%https://herbie.uwplse.org/demo/                         
%((L/sqrt(sqr(rho)+sqr(L-z)))-(z/sqrt(sqr(rho)+sqr(L-z)))+(z/sqrt(sqr(rho)+sqr(z))))
%(atan((L-z)./rho)+atan(z./rho))

% Reverse transform cyl to original Cart coord (trap centered ref)
%[Bout(:,1),Bout(:,2),Bout(:,3)]=pol2cart(theta-pi,B_perp,theta*0);
[Bout(:,1),Bout(:,2),Bout(:,3)]=pol2cart(theta-pi/2,B_perp,theta*0);

Bout=Bout*rev_rot_mat;


end