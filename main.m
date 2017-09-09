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
%   [x] trap generator
%   [x] trap (currents) --> trap freq, trap centre

verbose=1;

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
Iquad=1;
Ishunt=1.5;
Bbias=1e-4*[1,0,0];     % external bias field [T] (uniform assumption)

btrap=biquic_trap(Iquad,Ishunt,Bbias);

%%% Trap magnetic field calculation
[Bmag,Bxyz]=trap_eval(btrap,XYZ{:});

%% Trap visualisation - B-isosurfaces
hfig_btrap=plot_B_3d(btrap,Bmag,XYZ);

%% characterise trap
% evaluate trap center, 1D trap potential, trap frequencies
[f0,trap_cent,H]=trap_characterise(btrap,verbose);

%% End of code
t_end=toc(t_start);
disp('-----------------------------------------------');
fprintf('Total elapsed time (s): %7.1f\n',t_end);
disp('===================ALL TASKS COMPLETED===================');