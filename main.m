%% BiQUIC trap magnetic field calculator
% D K Shin
% refer to https://doi.org/10.1103/PhysRevA.35.1535 for expressions used

% physical constants
mu_0=4*pi*1e-7;     % vacuum permeability [Tm/A]

% TODO:
% [x] reverse transform mesghrid [xx,yy,zz] in trap coord to evaluate 
%   similarly ordered (TT,RR,ZZ) vectors for each unique coil
% [] evaluate B from each coil in trap coord
% [] sum B - total B field
% [] characterise potential
%   [] isosurfaces on B_tot

%% Config grid
% grid in trap centered ref frame
xyz0=cell(3,1);
xyz0{1}=linspace(-5e-3,5e-3,30);      % x-vect
xyz0{2}=linspace(-5e-3,5e-3,30);      % y-vect
xyz0{3}=linspace(0e-3,20e-3,30);      % z-vect

XYZ0=cell(3,1);
[XYZ0{1},XYZ0{2},XYZ0{3}]=meshgrid(xyz0{:});    % meshgrid

%% Scenario 1: a single coil
% config
R_coil=10e-3;    % coil radius [m]
I_coil=1;        % coil current [A]
ori_coil=[0,0,0];           % coil orientation - Euler angles
pos_coil=[0,0,0];           % coil centre position (x,y,z)

%%% evaluate the grid region in coil centered frame
% coord translation
XYZ0_coil=cell(3,1);
for ii=1:3
    XYZ0_coil{ii}=XYZ0{ii}+pos_coil(ii);     
end
% % cartesian to cylindrical
% [TT,RR,ZZ]=cart2pol(XYZ0_coil{:});
% 
% % evaluate the analytic magnetic field solution
% ell_k2=4*R_coil*RR./((R_coil+RR).^2+ZZ.^2);      % k^2 elliptic parameter
% [K,E]=ellipke(ell_k2);      % complete elliptic integrals of 1st,2nd orders
% 
% Bzz=(mu_0*I_coil./(2*pi)).*(1./sqrt((R_coil+RR).^2+ZZ.^2)).*(K+E.*(R_coil^2-RR.^2-ZZ.^2)./((R_coil-RR).^2+ZZ.^2));
% Brr=(mu_0*I_coil./(2*pi*RR)).*(ZZ./sqrt((R_coil+RR).^2+ZZ.^2)).*(-K+E.*(R_coil^2+RR.^2+ZZ.^2)./((R_coil-RR).^2+ZZ.^2));
% Bzz(~isfinite(Bzz))=NaN;    % Inf --> NaN
% Brr(~isfinite(Brr))=NaN;

% % Reverse transform cyl to original Cart coord (trap centered ref)
% [Bxx,Byy,Bzz]=pol2cart(TT,Brr,Bzz);

% Babs=sqrt(Bzz.^2+Brr.^2);     % absolute magnetic field strength [T]

[Bxx,Byy,Bzz]=Bfield_coil(R_coil,I_coil,XYZ0_coil{:});
Babs=sqrt(Bxx.^2+Byy.^2+Bzz.^2);     % absolute magnetic field strength [T]

%% Plot
% config
nBisosurf=6;

%%% Magnetic field
% quiver plot for B field
figure();
% quiver3(xx_tf,yy_tf,zz_tf,Bxx_tf,Byy_tf,Bzz_tf,...
%     'Color','k','LineWidth',1.5,'Visible','off');
quiver3(XYZ0{1},XYZ0{2},XYZ0{3},Bxx,Byy,Bzz,...
    'Color','k','LineWidth',1.5,'Visible','off');
hold on;

%%% B field magnitude - isosurfaces (isopotentials)
Bmax=max(Babs(isfinite(Babs)));
Bmin=min(Babs(isfinite(Babs)));
Bisoval=linspace(Bmin,Bmax,nBisosurf);
cc=viridis(nBisosurf);
p={};
for ii=1:nBisosurf
%     p{ii}=isosurface(xx_tf,yy_tf,zz_tf,Babs,Bisoval(ii));
    p{ii}=isosurface(XYZ0{1},XYZ0{2},XYZ0{3},Babs,Bisoval(ii));
    patch(p{ii},'FaceColor',cc(ii,:),'EdgeColor','none','FaceAlpha',0.15);
end
daspect([1,1,1]);
view(3);
axis tight;
camlight;
lighting gouraud;

xlabel('X');
ylabel('Y');
zlabel('Z');

%%% 2D B field magnitude contours (potential landscape)


%% functions
% useful rotator for meshgrid
function [x,y,z]=rotmesh(Mrot,x0,y0,z0)
nn=size(x0);    % size of the meshgrid
xyz=Mrot*[x0(:)';y0(:)';z0(:)'];
x=reshape(xyz(1,:),nn);
y=reshape(xyz(2,:),nn);
z=reshape(xyz(3,:),nn);
end

% B field calculator for single coil (located at origin, pointing Z-axis)
function [Bxx,Byy,Bzz]=Bfield_coil(R,I,x,y,z)
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
[Bxx,Byy,Bzz]=pol2cart(TT,Brr,Bzz);
end