%% BiQUIC trap magnetic field calculator
% D K Shin
% refer to https://doi.org/10.1103/PhysRevA.35.1535 for expressions used

clear all;

t_start=tic;

% TODO:
% [] Voltage to current conversion
% [] experiment params

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

%% Build BiQUIC trap
Iquad=1;
Ishunt=1.5;
Bbias=1e-4*[1,0,0];     % external bias field [T] (uniform assumption)

btrap=biquic_trap(Iquad,Ishunt,Bbias);  % build biquic

[Bmag,Bxyz]=trap_eval(btrap,XYZ{:});    % calculate magnetic field (macro summary)

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