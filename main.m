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
% [x] characterise potential
%   [x] isosurfaces on B magnitude
%   [x] trap frequency
% [] Voltage to current conversion
% [] experiment params
% [] package into a user friendly function
%   [] trap generator
%   [] trap (currents) --> trap freq, trap centre


%% Config grid
ngrid=50;          % 50 - med; 300 - very fine;

% grid in trap centered ref frame
xyz=cell(3,1);
%%% MACRO
xyz{1}=linspace(-20e-3,20e-3,ngrid);      % x-vect
xyz{2}=linspace(-20e-3,20e-3,ngrid);      % y-vect
xyz{3}=linspace(-20e-3,20e-3,ngrid);      % z-vect

% %%% MICRO
% xyz{1}=linspace(-1e-3,1e-3,ngrid);      % x-vect
% xyz{2}=linspace(-1e-3,1e-3,ngrid);      % y-vect
% xyz{3}=linspace(-1e-3,1e-3,ngrid);      % z-vect

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
Ishunt=1.5;

Rquad=Dquad/2;
Rshunt=Dshunt/2;

% Trap bias (nuller) [http://dx.doi.org/10.1063/1.2472600]
Bbias=1e-4*[1,0,0];     % external bias field [T] (uniform assumption)

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

%%% Trap magnetic field calculation
[Bmag,Bxyz]=trap_eval(btrap,XYZ{:});

%% Trap visualisation - B-isosurfaces
hfig_btrap=plot_B_3d(btrap,Bmag,XYZ);

%% Trap center
% trap center: point of minimum B magnitude

%%% Find trap center with current trap config
% find X to minimise $Bmag$ from trap_eval(btrap,X,0,0) - from symmetry
% initial guess param X~0.1 [mm]
% NOTE: x_cent solve in mm scale (function domain scaled to order of unity)
[x_cent,B_cent]=fminsearch(@(x) trap_eval(btrap,1e-3*x,0,0),0.1);   
trap_cent=[1e-3*x_cent,0,0];    % mm-->m evaluated trap centre [m]

% display evaluated trap centre
figure(hfig_btrap);
scatter3(1e3*trap_cent(1),1e3*trap_cent(2),1e3*trap_cent(3),...
    'Marker','o','MarkerFaceColor','r','MarkerEdgeColor','none',...
    'SizeData',30,'DisplayName','Trap centre');

%%% Get B field near trap centre
% config
ngrid_trap=50;
trap_lim=[-10e-6,10e-6; -10e-6,10e-6; -10e-6,10e-6];

% build trap scale grid (~10 um each dir)
xyz_trap=cell(3,1);
for ii=1:3
    xyz_trap{ii}=trap_cent(ii)+linspace(trap_lim(ii,1),trap_lim(ii,2),ngrid_trap);
end
XYZ_trap=cell(3,1);
[XYZ_trap{1},XYZ_trap{2},XYZ_trap{3}]=meshgrid(xyz_trap{:});    % meshgrid
% permute the 3D array so that indexing goes x-y-z
XYZ_trap=cellfun(@(YXZ) permute(YXZ,[2,1,3]),XYZ_trap,'UniformOutput',false);     

% calculate magnetic field
[Bmag_trap,Bxyz_trap]=trap_eval(btrap,XYZ_trap{:});

%%% 1D Bmag profile - X,Y,Z line profile
% indices to trap centre (in grid)
[~,I0]=min(abs(Bmag_trap(:)-B_cent));
II0=zeros(1,3);        % this is ordered in YXY
[II0(1),II0(2),II0(3)]=ind2sub(size(Bmag_trap),I0);     % index to trap cent

B_trap_1d=cell(3,1);     % 1D line-profile of magnetic field potential [T]
idxcirc=[1,2,3];
Bmagcirc=Bmag_trap;      % temporary copy
for ii=1:3
    B_trap_1d{ii}=Bmagcirc(:,II0(idxcirc(2)),II0(idxcirc(3)))';     % 1xNi array 
    idxcirc=circshift(idxcirc,[-1,0]);
    Bmagcirc=permute(Bmagcirc,[2,3,1]);     % dimension circular permutation
end
clearvars Bmagcirc;

%%% Visualise trap: 3D B-isosurfaces
hfig_btrap_cent_3d=plot_B_3d(btrap,Bmag_trap,XYZ_trap);

% display evaluated trap centre
figure(hfig_btrap_cent_3d);
scatter3(1e3*trap_cent(1),1e3*trap_cent(2),1e3*trap_cent(3),...
    'Marker','o','MarkerFaceColor','r','MarkerEdgeColor','none',...
    'SizeData',30,'DisplayName','Trap centre');

%%% Visualise B-1D: 1D B-profiles
hfig_bmag_1d=figure();
axisstr={'X','Y','Z'};
linestyle={'-','--',':'};
p=[];
xyz_trap0=cell(3,1);
for ii=1:3
    hold on;
    % shift coord to approximate trap center
    xyz_trap0{ii}=xyz_trap{ii}-trap_cent(ii);       % trap center ref'd disp [m]
    p(ii)=plot(1e3*xyz_trap0{ii},1e4*B_trap_1d{ii},...
        'LineStyle',linestyle{ii},'LineWidth',1.5,...
        'DisplayName',axisstr{ii});
end
box on;
xlabel('Displacement [mm]');
ylabel('$B$ [G]');
lgd=legend(p);
title(lgd,sprintf('(%0.2g,%0.2g,%0.2g) [mm]',1e3*trap_cent(:)));

%% Evaluate magnetic trap potential (1d: x-y-z)
% TODO - U_B± \propto mu_He * (±)Bmag (for weak(strong)-field seeking state)
const_prop=1;       % TODO
U_B_1d=cellfun(@(x) const_prop*x,B_trap_1d,'UniformOutput',false);


%% Trap frequency
% frequency: omega=sqrt(U''/m) [spatial 2nd order derivative]
% trap 1D: B_trap_1d - 3x1 cell array of 1D B magnitude profiles

% numerical second order differential - TODO check
d2U=cellfun(@(x)diff(x,2),U_B_1d,'UniformOutput',false);    % 2nd order DIFFERENCE in potential
dx=cellfun(@(x)diff(x),xyz_trap0,'UniformOutput',false);    % 1st ord DIFF in displacement

% evaluate trap frequency (harmonic)
m_He=1;     % TODO
f_trap=cellfun(@(D2Y,DX)(2*pi)*sqrt((D2Y./(DX(1:end-1).^2))/m_He),d2U,dx,'UniformOutput',false);     % trap frequency in [Hz]

f0=cellfun(@(x)mean(x),f_trap);     % mean trap frequency around trap centre [Hz]

%%% Visualise trap frequency - anharmonicity, etc
hfig_trap_freq=figure();
axisstr={'X','Y','Z'};
linestyle={'-','--',':'};
p=[];
for ii=1:3
    hold on;
    p(ii)=plot(1e3*xyz_trap0{ii}(2:end-1),f_trap{ii},...
        'LineStyle',linestyle{ii},'LineWidth',1.5,...
        'DisplayName',axisstr{ii});
end
box on;
xlabel('Displacement [mm]');
ylabel('$f$ [Hz]');
lgd=legend(p);


%% End of code
t_end=toc(t_start);
disp('-----------------------------------------------');
fprintf('Total elapsed time (s): %7.1f\n',t_end);
disp('===================ALL TASKS COMPLETED===================');