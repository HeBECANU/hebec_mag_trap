%% BiQUIC trap magnetic field calculator
% D K Shin
% refer to https://doi.org/10.1103/PhysRevA.35.1535 for expressions used

clear all;

t_start=tic;

% TODO:
% [x] reverse transform mesghrid [xx,yy,zz] in trap coord to evaluate 
%   similarly ordered (TT,RR,ZZ) vectors for each unique coil
% [x] evaluate B from each coil in trap coord
% [x] sum B - total B field
% [x] set up to biQUIC dimension - Ross
% [] Voltage to current conversion
% [] characterise potential
%   [] isosurfaces on B magnitude
%   [] trap frequency

%% Config grid
ngrid=50;          % 50 - med; 300 - very fine;

% grid in trap centered ref frame
xyz=cell(3,1);
% %%% MACRO
% xyz{1}=linspace(-20e-3,20e-3,ngrid);      % x-vect
% xyz{2}=linspace(-20e-3,20e-3,ngrid);      % y-vect
% xyz{3}=linspace(-20e-3,20e-3,ngrid);      % z-vect

%%% MICRO
xyz{1}=linspace(-1e-3,1e-3,ngrid);      % x-vect
xyz{2}=linspace(-1e-3,1e-3,ngrid);      % y-vect
xyz{3}=linspace(-1e-3,1e-3,ngrid);      % z-vect

XYZ=cell(3,1);
[XYZ{1},XYZ{2},XYZ{3}]=meshgrid(xyz{:});    % meshgrid
% permute the 3D array so that indexing goes x-y-z
XYZ=cellfun(@(YXZ) permute(YXZ,[2,1,3]),XYZ,'UniformOutput',false);     

%% Scenario 1: a single coil
% % config
% R_coil{1}=10e-3;    % coil radius [m]
% I_coil{1}=1;        % coil current [A]
% ori_coil{1}=[0,0,0];           % coil orientation - Euler angles
% pos_coil{1}=[0,0,0];           % coil centre position (x,y,z)
% 
% %%% Build trap as a struct array of components
% ncomps=numel(R_coil);
% objtype=cell(1,ncomps);
% objparam=cell(1,ncomps);
% objtype(:)={'coil'};    % all coils!
% for ii=1:ncomps
%     objparam{ii}={R_coil{ii},I_coil{ii},pos_coil{ii},ori_coil{ii}};
% end
% % create trap
% btrap=struct('type',objtype,'param',objparam);
% 
% %%% Trap magnetic field calculation
% [Bmag,Bxyz]=trap_eval(btrap,XYZ{:});

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
% %%% Build trap as a struct array of components
% ncomps=numel(R_coil);
% objtype=cell(1,ncomps);
% objparam=cell(1,ncomps);
% objtype(:)={'coil'};    % all coils!
% for ii=1:ncomps
%     objparam{ii}={R_coil{ii},I_coil{ii},pos_coil{ii},ori_coil{ii}};
% end
% btrap=struct('type',objtype,'param',objparam);
% 
% %%% Trap magnetic field calculation
% [Bmag,Bxyz]=trap_eval(btrap,XYZ{:});

%% Scenario 3: BiQUIC
%%% config
% Quad: 14mm Dia; 10 turns;
% Shunt (Ioffe): 14mm; 18 turns; 
% pitch=500um (wire dia); AH sep ~17mm; Q-Sh sep=18.5 mm; ~36 amp (max)
Dquad=14e-3;
Dshunt=14e-3;

Nturnshunt=18;
Nturnquad=10;
pitch_coil=0.5e-3;

disp_ah=17e-3;
disp_qs=18.5e-3;

Iquad=1;
Ishunt=0.2;

Rquad=Dquad/2;
Rshunt=Dshunt/2;

% Trap bias (nuller) [http://dx.doi.org/10.1063/1.2472600]
Bbias=1e-4*[0.01,0,0];     % external bias field (uniform assumption)

%%% Build trap
% Quadrupole - ref
quad_coil.type='coil';
quad_coil.param={Rquad,Iquad,[0,0,disp_ah/2],[0,0,0]};
% Shunt (Ioffe) - ref
shunt_coil.type='coil';
shunt_coil.param={Rshunt,Ishunt,[disp_qs,0,disp_ah/2],[0,0,0]};
% Bias (nuller)
bias.type='uniform';
bias.param={Bbias};

btrap=[];
for ii=1:Nturnquad
    quad_coil_temp=quad_coil;
    % shift by wire pitch
    quad_coil_temp.param{3}=quad_coil_temp.param{3}+(ii-1)*[0,0,pitch_coil];
    btrap=[btrap,quad_coil_temp];
    
    % anti-helmholtz pair - mirror symmetry around Z
    quad_coil_temp.param{2}=-1*quad_coil_temp.param{2};     % flip current dir
    quad_coil_temp.param{3}=[1,1,-1].*quad_coil_temp.param{3};
    btrap=[btrap,quad_coil_temp];
end
for ii=1:Nturnshunt
    shunt_coil_temp=shunt_coil;
    % shift by wire pitch
    shunt_coil_temp.param{3}=shunt_coil_temp.param{3}+(ii-1)*[0,0,pitch_coil];
    btrap=[btrap,shunt_coil_temp];
    
    % anti-helmholtz pair - mirror symmetry around Z
    shunt_coil_temp.param{2}=-1*shunt_coil_temp.param{2};   % flip current dir
    shunt_coil_temp.param{3}=[1,1,-1].*shunt_coil_temp.param{3};
    btrap=[btrap,shunt_coil_temp];
end
btrap=[btrap,bias];

%%% Find trap center with current trap config
% find X to minimise $Bmag$ from trap_eval(btrap,X,0,0) - from symmetry
% initial guess param ~1e-6

%%% Trap magnetic field calculation
[Bmag,Bxyz]=trap_eval(btrap,XYZ{:});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EDIT - bias incorporated in trap definition, vector field packaged in a
% single cell-array output (see Bfield_coil.m and trap_eval.m)
% 
% % apply bias field (nuller; spatially uniform)
% % http://dx.doi.org/10.1063/1.2472600
% Bxx=Bxx+Bbias(1);
% Byy=Byy+Bbias(2);
% Bzz=Bzz+Bbias(3);
% 
% Bmag=sqrt(Bxx.^2+Byy.^2+Bzz.^2);     % absolute magnetic field strength [T]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Plot
% config
nBisosurf=5;

%%% Magnetic field
% quiver plot for B field
hfig_btrap=figure();
% quiver3(1e3*XYZ{1},1e3*XYZ{2},1e3*XYZ{3},Bxyz{:},...
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
    p{ii}=isosurface(1e3*XYZ{1},1e3*XYZ{2},1e3*XYZ{3},Bmag,Bisoval(ii));
    pp(ii)=patch(p{ii},'FaceColor',cc(ii,:),'EdgeColor','none','FaceAlpha',0.15,...
        'DisplayName',sprintf('%0.1g',1e4*Bisoval(ii)));
end
box on;
daspect([1,1,1]);
view(3);
xlim(1e3*[min(xyz{1}),max(xyz{1})]);
ylim(1e3*[min(xyz{2}),max(xyz{2})]);
zlim(1e3*[min(xyz{3}),max(xyz{3})]);
camlight;
lighting gouraud;

xlabel('X [mm]');
ylabel('Y [mm]');
zlabel('Z [mm]');

%%% Coils
% NOTE: coil widths are not to true scale!
% unit ring
nnring=100;
phi=linspace(0,2*pi,nnring);
[x_ring,y_ring]=pol2cart(phi,1);
z_ring=zeros(1,nnring);
figure(hfig_btrap);
for ii=1:numel(btrap)
    % only plot coils
    if isequal(btrap(ii).type,'coil')
        % transform from unit ring
        R_coil_this=btrap(ii).param{1};
        pos_coil_this=btrap(ii).param{3};
        xthis=R_coil_this*x_ring+pos_coil_this(1);
        ythis=R_coil_this*y_ring+pos_coil_this(2);
        zthis=R_coil_this*z_ring+pos_coil_this(3);
        
        % draw this coil
        hold on;
        plot3(1e3*xthis,1e3*ythis,1e3*zthis,...
            'Color','k','LineWidth',2);
    end
end

% legend
lgd=legend(pp(:));
title(lgd,'$B$-isosurface (G)');

%%% 2D B field magnitude contours (potential landscape)

%% Characterise trap: trap center and frequency
% trap center: point of minimum B magnitude
% frequency: omega=sqrt(V''/m) [spatial 2nd order derivative]

%%% approximate trap centre (evaluated grid)
[B0_approx,I0_approx]=min(Bmag(:));
xyz0_approx=cellfun(@(x) x(I0_approx),XYZ);
II0_approx=zeros(1,3);        % this is ordered in YXY
[II0_approx(1),II0_approx(2),II0_approx(3)]=ind2sub(size(Bmag),I0_approx);

%%% B profile in X,Y,Z line profile
B_1d=cell(3,1);     % 1D line-profile of magnetic field potential [T]
idxcirc=[1,2,3];
Bmagcirc=Bmag;      % temporary copy
for ii=1:3
    B_1d{ii}=Bmagcirc(:,II0_approx(idxcirc(2)),II0_approx(idxcirc(3)));
    idxcirc=circshift(idxcirc,-1);
    Bmagcirc=permute(Bmagcirc,[2,3,1]);     % dimension circular permutation 
end
clearvars Bmagcirc;

% plot
hfig_bmag_1d=figure();
axisstr={'X','Y','Z'};
linestyle={'-','--',':'};
p=[];
for ii=1:3
    hold on;
    % shift coord to approximate trap center
    p(ii)=plot(1e3*(xyz{ii}-xyz0_approx(ii)),1e4*B_1d{ii},...
        'LineStyle',linestyle{ii},'LineWidth',1.5,...
        'DisplayName',axisstr{ii});
end
box on;
xlabel('Displacement [mm]');
ylabel('$B$ [G]');
lgd=legend(p);
title(lgd,sprintf('(%0.2g,%0.2g,%0.2g) [mm]',1e3*xyz0_approx(:)));

%%% approximate trap frequency (evaluated grid)


%% End of code
t_end=toc(t_start);
disp('-----------------------------------------------');
fprintf('Total elapsed time (s): %7.1f\n',t_end);
disp('===================ALL TASKS COMPLETED===================');