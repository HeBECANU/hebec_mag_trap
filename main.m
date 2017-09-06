%% BiQUIC trap magnetic field calculator
% D K Shin
% refer to https://doi.org/10.1103/PhysRevA.35.1535 for expressions used

t_start=tic;

% TODO:
% [x] reverse transform mesghrid [xx,yy,zz] in trap coord to evaluate 
%   similarly ordered (TT,RR,ZZ) vectors for each unique coil
% [x] evaluate B from each coil in trap coord
% [x] sum B - total B field
% [] set up to biQUIC dimension - Ross
% [] Voltage to current conversion
% [] characterise potential
%   [] isosurfaces on B magnitude
%   [] trap frequency

%% Config grid
ngrid=100;          % 50 - med; 300 - very fine;

% grid in trap centered ref frame
xyz0=cell(3,1);
%%% MACRO
xyz0{1}=linspace(-20e-3,20e-3,ngrid);      % x-vect
xyz0{2}=linspace(-20e-3,20e-3,ngrid);      % y-vect
xyz0{3}=linspace(-20e-3,20e-3,ngrid);      % z-vect

% %%% MICRO
% xyz0{1}=linspace(-1e-3,1e-3,ngrid);      % x-vect
% xyz0{2}=linspace(-1e-3,1e-3,ngrid);      % y-vect
% xyz0{3}=linspace(-1e-3,1e-3,ngrid);      % z-vect

XYZ0=cell(3,1);
[XYZ0{1},XYZ0{2},XYZ0{3}]=meshgrid(xyz0{:});    % meshgrid

%% Scenario 1: a single coil
% % config
% R_coil=10e-3;    % coil radius [m]
% I_coil=1;        % coil current [A]
% ori_coil=[0,0,0];           % coil orientation - Euler angles
% pos_coil=[0,0,0];           % coil centre position (x,y,z)
% 
% %%% evaluate the grid region in coil centered frame
% % coord translation
% XYZ0_coil=cell(3,1);
% for ii=1:3
%     XYZ0_coil{ii}=XYZ0{ii}-pos_coil(ii);     
% end
% 
% [Bxx,Byy,Bzz]=Bfield_coil(R_coil,I_coil,XYZ0_coil{:});
% Bmag=sqrt(Bxx.^2+Byy.^2+Bzz.^2);     % absolute magnetic field strength [T]

%% Scenario 2: anti-Helmholtz with 2 coils
% %%% config
% % coil 1
% R_coil{1}=10e-3;    % coil radius [m]
% I_coil{1}=-1;        % coil current [A]
% ori_coil{1}=[0,0,0];           % coil orientation - Euler angles
% pos_coil{1}=[0,0,-4e-3];           % coil centre position (x,y,z)
% 
% % coil 2
% R_coil{2}=10e-3;    % coil radius [m]
% I_coil{2}=1;        % coil current [A]
% ori_coil{2}=[0,0,0];           % coil orientation - Euler angles
% pos_coil{2}=[0,0,4e-3];           % coil centre position (x,y,z)
% 
% %%% evaluate the grid region in coil centered frame
% for ii=1:length(R_coil)
%     R_coil_this=R_coil{ii};
%     I_coil_this=I_coil{ii};
%     pos_coil_this=pos_coil{ii};
%     
%     % coord translation
%     XYZ0_coil_this=cell(3,1);
%     for jj=1:3
%         XYZ0_coil_this{jj}=XYZ0{jj}-pos_coil_this(jj);
%     end
%     
%     Bxyz_this=cell(3,1);
%     [Bxyz_this{1},Bxyz_this{2},Bxyz_this{3}]=Bfield_coil(R_coil_this,I_coil_this,XYZ0_coil_this{:});
%     
%     if ii==1
%         Bxyz=Bxyz_this;
%     else
%         Bxyz=cellfun(@(x,y)x+y,Bxyz,Bxyz_this,'UniformOutput',false);
%     end
% end
% Bxx=Bxyz{1};
% Byy=Bxyz{2};
% Bzz=Bxyz{3};
% 
% Bmag=sqrt(Bxx.^2+Byy.^2+Bzz.^2);     % absolute magnetic field strength [T]

%% Scenario 3: BiQUIC
%%% config
Rquad=5e-3;
Iquad=1;
Rshunt=3e-3;
Ishunt=1;

% QUAD coil 1
R_coil{1}=Rquad;    % coil radius [m]
I_coil{1}=-Iquad;        % coil current [A]
ori_coil{1}=[0,0,0];           % coil orientation - Euler angles
pos_coil{1}=[0,0,-4e-3];           % coil centre position (x,y,z)

% Quad coil 2
R_coil{2}=Rquad;    % coil radius [m]
I_coil{2}=Iquad;        % coil current [A]
ori_coil{2}=[0,0,0];           % coil orientation - Euler angles
pos_coil{2}=[0,0,4e-3];           % coil centre position (x,y,z)

% Shunt
% Shunt coil 1
R_coil{3}=Rshunt;    % coil radius [m]
I_coil{3}=-Ishunt;        % coil current [A]
ori_coil{3}=[0,0,0];           % coil orientation - Euler angles
pos_coil{3}=[10e-3,0,-4e-3];           % coil centre position (x,y,z)

% Shunt coil 2
R_coil{4}=Rshunt;    % coil radius [m]
I_coil{4}=Ishunt;        % coil current [A]
ori_coil{4}=[0,0,0];           % coil orientation - Euler angles
pos_coil{4}=[10e-3,0,4e-3];           % coil centre position (x,y,z)

%%% evaluate the grid region in coil centered frame
for ii=1:length(R_coil)
    R_coil_this=R_coil{ii};
    I_coil_this=I_coil{ii};
    pos_coil_this=pos_coil{ii};
    
    % coord translation
    XYZ0_coil_this=cell(3,1);
    for jj=1:3
        XYZ0_coil_this{jj}=XYZ0{jj}-pos_coil_this(jj);
    end
    
    Bxyz_this=cell(3,1);
    [Bxyz_this{1},Bxyz_this{2},Bxyz_this{3}]=Bfield_coil(R_coil_this,I_coil_this,XYZ0_coil_this{:});
    
    if ii==1
        Bxyz=Bxyz_this;
    else
        Bxyz=cellfun(@(x,y)x+y,Bxyz,Bxyz_this,'UniformOutput',false);
    end
end
Bxx=Bxyz{1};
Byy=Bxyz{2};
Bzz=Bxyz{3};

Bmag=sqrt(Bxx.^2+Byy.^2+Bzz.^2);     % absolute magnetic field strength [T]

%% Plot
% config
nBisosurf=4;

%%% Magnetic field
% quiver plot for B field
hfig_btrap=figure();
% quiver3(1e3*XYZ0{1},1e3*XYZ0{2},1e3*XYZ0{3},Bxx,Byy,Bzz,...
%     'Color','k','LineWidth',1.5,'Visible','off');
hold on;

%%% B field magnitude - isosurfaces (isopotentials)
Bmax=max(Bmag(isfinite(Bmag)));
Bmin=min(Bmag(isfinite(Bmag)));
% Bisoval=linspace(Bmin,Bmax,nBisosurf+2);
Bisoval=logspace(log10(Bmin),log10(Bmax),nBisosurf+2);
Bisoval=Bisoval(2:end-1);   % cull the min and max
cc=viridis(nBisosurf);
p={};
pp=[];
for ii=1:nBisosurf
    p{ii}=isosurface(1e3*XYZ0{1},1e3*XYZ0{2},1e3*XYZ0{3},Bmag,Bisoval(ii));
    pp(ii)=patch(p{ii},'FaceColor',cc(ii,:),'EdgeColor','none','FaceAlpha',0.15,...
        'DisplayName',sprintf('%0.1g',1e4*Bisoval(ii)));
end
box on;
daspect([1,1,1]);
view(3);
xlim(1e3*[min(xyz0{1}),max(xyz0{1})]);
ylim(1e3*[min(xyz0{2}),max(xyz0{2})]);
zlim(1e3*[min(xyz0{3}),max(xyz0{3})]);
camlight;
lighting gouraud;

lgd=legend(pp(:));
title(lgd,'$B$-isosurface (G)');

xlabel('X [mm]');
ylabel('Y [mm]');
zlabel('Z [mm]');

%%% 2D B field magnitude contours (potential landscape)


%% End of code
t_end=toc(t_start);
disp('-----------------------------------------------');
fprintf('Total elapsed time (s): %7.1f\n',t_end);
disp('===================ALL TASKS COMPLETED===================');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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