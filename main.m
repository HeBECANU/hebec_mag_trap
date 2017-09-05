%% BiQUIC trap magnetic field calculator
% D K Shin
% refer to https://doi.org/10.1103/PhysRevA.35.1535 for expressions used

% physical constants
mu_0=4*pi*1e-7;     % vacuum permeability [Tm/A]

% TODO:
% [] reverse transform mesghrid [xx,yy,zz] in trap coord to evaluate 
%   similarly ordered (TT,RR,ZZ) vectors for each unique coil
% [] evaluate B from each coil in trap coord
% [] sum B - total B field
% [] characterise potential
%   [] isosurfaces on B_tot

%% Scenario 1: a single coil
% config
R=10e-3;    % coil radius [m]
I=1;        % coil current [A]
ori=[0,0,0];           % coil orientation - Euler angles
pos=[0,0,-10e-3];      % coil centre position (x,y,z)


% region to evaluate
nzz=20;
nrr=20;
ntheta=20;
zz=linspace(-20e-3,20e-3,nzz);    % axial [m]
rr=linspace(0,20e-3,nrr);        % radial (cylinder) [m]
theta=linspace(0,2*pi,ntheta);  % azim angle [rad]
[TT,RR,ZZ]=meshgrid(theta,rr,zz);

% evaluate the analytic magnetic field solution
ell_k2=4*R*RR./((R+RR).^2+ZZ.^2);      % k^2 elliptic parameter
[K,E]=ellipke(ell_k2);      % complete elliptic integrals of 1st,2nd orders

Bzz=(mu_0*I./(2*pi)).*(1./sqrt((R+RR).^2+ZZ.^2)).*(K+E.*(R^2-RR.^2-ZZ.^2)./((R-RR).^2+ZZ.^2));
Brr=(mu_0*I./(2*pi*RR)).*(ZZ./sqrt((R+RR).^2+ZZ.^2)).*(-K+E.*(R^2+RR.^2+ZZ.^2)./((R-RR).^2+ZZ.^2));
Bzz(~isfinite(Bzz))=NaN;    % Inf --> NaN
Brr(~isfinite(Brr))=NaN;

Babs=sqrt(Bzz.^2+Brr.^2);     % absolute magnetic field strength [T]

% evaluate in Cartesian coordinate system
[xx,yy,zz]=pol2cart(TT,RR,ZZ);
[Bxx,Byy,Bzz]=pol2cart(TT,Brr,Bzz);

% coil orientation + position transform
Mrot=euler2rotm(ori);
[xx_tf,yy_tf,zz_tf]=rotmesh(Mrot,xx,yy,zz);
[Bxx_tf,Byy_tf,Bzz_tf]=rotmesh(Mrot,Bxx,Byy,Bzz);

xx_tf=xx_tf+pos(1);
yy_tf=yy_tf+pos(2);
zz_tf=zz_tf+pos(3);

% plot the magnetic field
% quiver plot for B field
figure();
quiver3(xx_tf,yy_tf,zz_tf,Bxx_tf,Byy_tf,Bzz_tf,...
    'Color','k','LineWidth',1.5,'Visible','off');
hold on;

nBisosurf=6;
Bmax=max(Babs(isfinite(Babs)));
Bmin=min(Babs(isfinite(Babs)));
Bisoval=linspace(Bmin,Bmax,nBisosurf);
cc=viridis(nBisosurf);
p={};
for ii=1:nBisosurf
    p{ii}=isosurface(xx_tf,yy_tf,zz_tf,Babs,Bisoval(ii));
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

% Babs contours (potential landscape)


%% functions
% useful rotator for meshgrid
function [x,y,z]=rotmesh(Mrot,x0,y0,z0)
nn=size(x0);    % size of the meshgrid
xyz=Mrot*[x0(:)';y0(:)';z0(:)'];
x=reshape(xyz(1,:),nn);
y=reshape(xyz(2,:),nn);
z=reshape(xyz(3,:),nn);
end