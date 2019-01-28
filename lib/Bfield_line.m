function Bout=Bfield_line(R,I,rot_vec,xyz)
% Bout = Bfield_coil(R, I, x, y, z)
% B field calculator for line of current (starting at origin, pointing Z-axis)
% https://www.miniphysics.com/uy1-magnetic-field-of-a-straight-current-carrying-conductor.html
% Bout is a 3x1 cell-array: {Bx,By,Bz} defined at points x,y,z

%hack to get to work in parfor
mu_0=1.2566370614*10^-6;

rot_mat=rotationVectorToMatrix(rot_vec);
xyz=xyz*rot_mat;
% cartesian to cylindrical
% theta-rad-z
[TT,RR,ZZ]=cart2pol(xyz(:,1),xyz(:,2),xyz(:,3));



error('not yet implemented')
Bzz=0;
Brr=(mu_0*I./(2*pi*RR)).*(ZZ./sqrt((R+RR).^2+ZZ.^2)).*(-K+E.*(R^2+RR.^2+ZZ.^2)./((R-RR).^2+ZZ.^2));
Bzz(~isfinite(Bzz))=NaN;    % Inf --> NaN
Brr(~isfinite(Brr))=NaN;

% Reverse transform cyl to original Cart coord (trap centered ref)

[Bout(:,1),Bout(:,2),Bout(:,3)]=pol2cart(TT,Brr,Bzz);
end