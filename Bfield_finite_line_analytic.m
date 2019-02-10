function Bout=Bfield_finite_line_analytic(tlim,curr,rot_vec,xyz)
% Bout = Bfield_finite_line_analytic
% B field calculator for single line of current (located at origin, pointing Z-axis)
% Bout is a 3x1 cell-array: {Bx,By,Bz} defined at points x,y,z
% 20190210 validated to give the same answer as Bfield_path_numeric 
% TO DO

% physical constants
%global const
%mu_0=const.mu0;  
%hack to get to work in parfor
mu_0=1.2566370614*10^-6;

rot_mat=rotationVectorToMatrix(rot_vec);
rev_rot_mat=rotationVectorToMatrix(-rot_vec); %reverse rotation vector
xyz=xyz*rot_mat;

%for derivation of the expression see mathematica notebook mag_field_of_line_and_helix
%also see %http://web.mit.edu/viz/EM/visualizations/coursenotes/modules/guide09.pdf
%x=(py ((pz-tmax)/Sqrt[px^2+py^2+(pz-tmax)^2]+(-pz+tmin)/Sqrt[px^2+py^2+(pz-tmin)^2]))/(px^2+py^2)
%y=(px ((-pz+tmax)/Sqrt[px^2+py^2+(pz-tmax)^2]+(pz-tmin)/Sqrt[px^2+py^2+(pz-tmin)^2]))/(px^2+py^2)
%z=0

Bout=[(xyz(:,2)./(xyz(:,1).^2 + xyz(:,2).^2)).*(...
                (xyz(:,3)-tlim(2))./sqrt(xyz(:,1).^2 +xyz(:,2).^2+(xyz(:,3)-tlim(2)).^2) + ...
                (-xyz(:,3)+tlim(1))./sqrt(xyz(:,1).^2+xyz(:,2).^2+(xyz(:,3)-tlim(1)).^2))...
    ,(xyz(:,1)./(xyz(:,1).^2 + xyz(:,2).^2)).*(...
                (-xyz(:,3)+tlim(2))./sqrt(xyz(:,1).^2 +xyz(:,2).^2+(xyz(:,3)-tlim(2)).^2) + ...
                (xyz(:,3)-tlim(1))./sqrt(xyz(:,1).^2+xyz(:,2).^2+(xyz(:,3)-tlim(1)).^2))...
    ,xyz(:,3)*0];

Bout=((mu_0*curr)/(4*pi)).*Bout;

Bout=Bout*rev_rot_mat;
end

