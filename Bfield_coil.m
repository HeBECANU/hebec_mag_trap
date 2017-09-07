function Bout=Bfield_coil(R,I,x,y,z)
% Bout = Bfield_coil(R, I, x, y, z)
% B field calculator for single coil (located at origin, pointing Z-axis)
% Bout is a 3x1 cell-array: {Bx,By,Bz} defined at points x,y,z


% physical constants
mu_0=4*pi*1e-7;     % vacuum permeability [Tm/A]

% cartesian to cylindrical
% theta-rad-z
[TT,RR,ZZ]=cart2pol(x,y,z);

% evaluate the analytic magnetic field solution
% https://doi.org/10.1103/PhysRevA.35.1535
ell_k2=4*R*RR./((R+RR).^2+ZZ.^2);      % k^2 elliptic parameter
[K,E]=ellipke(ell_k2);      % complete elliptic integrals of 1st,2nd orders

Bzz=(mu_0*I./(2*pi)).*(1./sqrt((R+RR).^2+ZZ.^2)).*(K+E.*(R^2-RR.^2-ZZ.^2)./((R-RR).^2+ZZ.^2));
Brr=(mu_0*I./(2*pi*RR)).*(ZZ./sqrt((R+RR).^2+ZZ.^2)).*(-K+E.*(R^2+RR.^2+ZZ.^2)./((R-RR).^2+ZZ.^2));
Bzz(~isfinite(Bzz))=NaN;    % Inf --> NaN
Brr(~isfinite(Brr))=NaN;

% Reverse transform cyl to original Cart coord (trap centered ref)
Bout=cell(3,1);
[Bout{:}]=pol2cart(TT,Brr,Bzz);
end